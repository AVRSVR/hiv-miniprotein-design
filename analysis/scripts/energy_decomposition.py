"""
Is D5's interface electrostatics-dominated compared with real gp120 complexes?

D5 is an E/K-rich helical hairpin with a GRAVY of -1.50 and a single salt bridge,
which is the fingerprint of a binder that sticks to charged surfaces rather than
recognising one. The way to test that from the trajectories we already have is to
split the MM/GBSA binding energy into its van der Waals part and everything
electrostatic (Coulomb plus polar solvation), and compare the split against the
experimental CD4-gp120 and b12-gp120 interfaces scored the same way.

The split is obtained by scoring each frame twice: once normally, and once with
every partial charge zeroed in both the nonbonded and the GB force, which leaves
the Lennard-Jones term alone. The difference is the electrostatic contribution.

    dG_total = dG_vdW + dG_elec+solv

A design whose interface is far more electrostatics-weighted than the natural
complexes is a specificity risk; one in the same range is not.
"""
import os, sys, json, warnings
import numpy as np
warnings.filterwarnings("ignore")

import mdtraj as mdt
from openmm.app import PDBFile
from pdbfixer import PDBFixer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmgbsa import build_topology, make_context, subset, energy, SPLIT, MD, OUT

TRAJ = {"S310": "D1 (CD4-mimetic control)", "Chai1_Best": "D2 (de novo bispecific)",
        "Cyclic_Ring": "D3 (CD4-mimetic cyclic)", "model_4_chai": "D2' (alternative pose)",
        "Latest_Best": "D5 (de novo, coreceptor site)"}
# experimental references, scored as single minimized structures
REF = {"1GC1": ("CD4-gp120 (experimental)", "G", "C"),
       "2NY7": ("b12-gp120 (experimental)", "G", "H")}
PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def zero_charges(sim):
    """Strip every partial charge, leaving Lennard-Jones only."""
    sysm = sim.context.getSystem()
    for f in sysm.getForces():
        n = f.__class__.__name__
        if n == "NonbondedForce":
            for i in range(f.getNumParticles()):
                q, s, e = f.getParticleParameters(i)
                f.setParticleParameters(i, 0.0, s, e)
            for i in range(f.getNumExceptions()):
                p1, p2, qq, s, e = f.getExceptionParameters(i)
                f.setExceptionParameters(i, p1, p2, 0.0, s, e)
            f.updateParametersInContext(sim.context)
        elif n == "CustomGBForce":
            for i in range(f.getNumParticles()):
                p = list(f.getParticleParameters(i))
                p[0] = 0.0
                f.setParticleParameters(i, p)
            f.updateParametersInContext(sim.context)
    return sim


def split_energies(top, pos, frames):
    """(total, vdW) binding energies averaged over frames, in kJ/mol."""
    atoms = list(top.atoms())
    rec = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "A"])
    lig = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "B"])
    tr, pr = subset(top, pos, "A")
    tl, pl = subset(top, pos, "B")
    full = [make_context(top, pos), make_context(tr, pr), make_context(tl, pl)]
    lj = [zero_charges(make_context(top, pos)),
          zero_charges(make_context(tr, pr)),
          zero_charges(make_context(tl, pl))]
    tot, vdw = [], []
    for x in frames:
        tot.append(energy(full[0], x) - energy(full[1], x[rec]) - energy(full[2], x[lig]))
        vdw.append(energy(lj[0], x) - energy(lj[1], x[rec]) - energy(lj[2], x[lig]))
    return np.array(tot), np.array(vdw)


def row(label, tot, vdw):
    t, v = float(tot.mean()), float(vdw.mean())
    e = t - v
    return {"system": label,
            "dG_total_kcal": round(t / 4.184, 1),
            "dG_vdW_kcal": round(v / 4.184, 1),
            "dG_elec_solv_kcal": round(e / 4.184, 1),
            "elec_fraction": round(e / t, 3) if t else None,
            "n_frames": len(tot)}


rows = []
for tag, label in TRAJ.items():
    sp = os.path.join(SPLIT, f"{tag}_split.pdb")
    dcd = os.path.join(MD, f"{tag}_traj.dcd")
    if not (os.path.exists(sp) and os.path.exists(dcd)):
        continue
    top, pos = build_topology(sp)
    ref = os.path.join(MD, f"{tag}_decomp_top.pdb")
    with open(ref, "w") as fh:
        PDBFile.writeFile(top, pos, fh)
    tr = mdt.load(dcd, top=ref)
    frames = tr.xyz[len(tr) // 2:][::3]          # final half, every third frame
    tot, vdw = split_energies(top, pos, frames)
    rows.append(row(label, tot, vdw))
    print(f"[ok] {label:<32} total {rows[-1]['dG_total_kcal']:>8.1f} "
          f"vdW {rows[-1]['dG_vdW_kcal']:>8.1f} elec {rows[-1]['dG_elec_solv_kcal']:>8.1f} "
          f"| elec frac {rows[-1]['elec_fraction']}", flush=True)

for pdb, (label, rc, lc) in REF.items():
    src = os.path.join(PA, "01_structures", f"{pdb}.pdb")
    if not os.path.exists(src):
        continue
    fx = PDBFixer(filename=src)
    keep = {rc, lc}
    for ch in list(fx.topology.chains()):
        pass
    fx.removeChains([i for i, ch in enumerate(fx.topology.chains()) if ch.id not in keep])
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)
    # relabel to the A/B convention the splitter uses
    for ch in fx.topology.chains():
        ch.id = "A" if ch.id == rc else "B"
    sim = make_context(fx.topology, fx.positions)
    sim.minimizeEnergy(maxIterations=800)
    st = sim.context.getState(getPositions=True)
    import openmm.unit as u
    xyz = np.array(st.getPositions().value_in_unit(u.nanometer))
    tot, vdw = split_energies(fx.topology, st.getPositions(), np.array([xyz]))
    rows.append(row(label, tot, vdw))
    print(f"[ok] {label:<32} total {rows[-1]['dG_total_kcal']:>8.1f} "
          f"vdW {rows[-1]['dG_vdW_kcal']:>8.1f} elec {rows[-1]['dG_elec_solv_kcal']:>8.1f} "
          f"| elec frac {rows[-1]['elec_fraction']}", flush=True)

with open(os.path.join(OUT, "energy_decomposition.json"), "w") as fh:
    json.dump(rows, fh, indent=2)
print("DECOMPOSITION DONE", flush=True)
