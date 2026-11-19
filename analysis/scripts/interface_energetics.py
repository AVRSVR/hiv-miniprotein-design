"""
Per-residue interface energetics across the MD ensemble.

Decomposes the binder-receptor interaction into per-residue contributions, so the
interface can be described by which residues actually carry it rather than by bulk
totals. For every binder residue and every sampled frame we sum its Lennard-Jones
and electrostatic interaction with all receptor atoms, using the force field's own
parameters, and average over the trajectory.

Electrostatics use a distance-dependent dielectric, eps(r) = 4r, the standard
screening approximation for interface decomposition in implicit solvent; absolute
values therefore carry the usual caveats, but the *ranking* of residues - which is
what a hotspot analysis is for - is robust to the choice.

Also reports per-residue buried surface area, averaged over the same frames. A
residue that is both energetically dominant and extensively buried is a hotspot:
the position a redesign must preserve and an escape mutation would target.

Outputs 05_tables/interface_energetics.json
"""
import os, sys, json, glob, warnings
import numpy as np
warnings.filterwarnings("ignore")

import mdtraj as md
from openmm.app import ForceField, HBonds, NoCutoff, PDBFile
from openmm.unit import amu, nanometer
from pdbfixer import PDBFixer

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SPLIT = os.path.join(PROJECT, "paper_analysis", "01_structures")
MD = os.path.join(PROJECT, "paper_analysis", "03_md")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")

COULOMB = 138.935458          # kJ nm / (mol e^2)
CUTOFF = 1.2                  # nm
STRIDE = 5                    # every 5th frame (50 ps)

LABEL = {"S310": "D1", "Chai1_Best": "D2", "Cyclic_Ring": "D3",
         "model_4_chai": "D4", "Latest_Best": "D5"}


def build(split_pdb):
    fx = PDBFixer(filename=split_pdb)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)
    return fx.topology, fx.positions


def nb_params(topology):
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(topology, nonbondedMethod=NoCutoff,
                             constraints=HBonds, hydrogenMass=1.5 * amu)
    nbf = None
    for f in system.getForces():
        if f.__class__.__name__ == "NonbondedForce":
            nbf = f
            break
    if nbf is None:
        raise SystemExit("no NonbondedForce")
    n = system.getNumParticles()
    q = np.zeros(n); sig = np.zeros(n); eps = np.zeros(n)
    for i in range(n):
        c, s, e = nbf.getParticleParameters(i)
        q[i] = c.value_in_unit_system(__import__("openmm").unit.md_unit_system)
        sig[i] = s.value_in_unit(nanometer)
        eps[i] = e.value_in_unit_system(__import__("openmm").unit.md_unit_system)
    return q, sig, eps


def run(tag):
    split_pdb = os.path.join(SPLIT, f"{tag}_split.pdb")
    dcd = os.path.join(MD, f"{tag}_traj.dcd")
    if not (os.path.exists(split_pdb) and os.path.exists(dcd)):
        return None

    top, pos = build(split_pdb)
    ref = os.path.join(MD, f"{tag}_mmgbsa_top.pdb")
    if not os.path.exists(ref):
        with open(ref, "w") as fh:
            PDBFile.writeFile(top, pos, fh)

    q, sig, eps = nb_params(top)
    atoms = list(top.atoms())
    rec = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "A"])
    lig = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "B"])
    if rec.size == 0 or lig.size == 0:
        return None

    traj = md.load(dcd, top=ref)
    frames = range(0, traj.n_frames, STRIDE)

    lig_res = sorted({atoms[i].residue.index for i in lig})
    res_atoms = {r: np.array([i for i in lig if atoms[i].residue.index == r])
                 for r in lig_res}
    res_name = {r: (atoms[res_atoms[r][0]].residue.name,
                    atoms[res_atoms[r][0]].residue.id) for r in lig_res}

    acc = {r: [] for r in lig_res}
    for f in frames:
        x = traj.xyz[f]
        xr = x[rec]
        for r in lig_res:
            ai = res_atoms[r]
            d = np.linalg.norm(x[ai][:, None, :] - xr[None, :, :], axis=-1)
            m = d < CUTOFF
            if not m.any():
                acc[r].append(0.0)
                continue
            ii, jj = np.where(m)
            rr = d[ii, jj]
            qi = q[ai][ii]; qj = q[rec][jj]
            si = (sig[ai][ii] + sig[rec][jj]) / 2.0
            ei = np.sqrt(eps[ai][ii] * eps[rec][jj])
            sr6 = (si / rr) ** 6
            lj = 4.0 * ei * (sr6 ** 2 - sr6)
            # distance-dependent dielectric eps(r) = 4r
            coul = COULOMB * qi * qj / (4.0 * rr * rr)
            acc[r].append(float((lj + coul).sum()))

    # per-residue buried surface area over the same frames
    sub = traj[list(frames)]
    t = sub.topology
    rec_sel = t.select("chainid 0")
    lig_sel = t.select("chainid 1")
    sasa_c = md.shrake_rupley(sub, mode="residue")
    sasa_l = md.shrake_rupley(sub.atom_slice(lig_sel), mode="residue")
    n_rec_res = len({t.atom(int(i)).residue.index for i in rec_sel})

    rows = []
    for k, r in enumerate(lig_res):
        e = np.array(acc[r])
        bsa = float((sasa_l[:, k].mean() - sasa_c[:, n_rec_res + k].mean()) * 100.0)
        nm, rid = res_name[r]
        rows.append({"resname": nm, "resid": int(rid),
                     "energy_kJ_mol": round(float(e.mean()), 2),
                     "energy_sd": round(float(e.std()), 2),
                     "buried_SASA_A2": round(bsa, 1)})
    rows.sort(key=lambda z: z["energy_kJ_mol"])
    total = sum(z["energy_kJ_mol"] for z in rows)
    return {"design": tag, "label": LABEL.get(tag, tag),
            "n_frames": len(list(frames)), "total_kJ_mol": round(total, 1),
            "residues": rows}


if __name__ == "__main__":
    tags = sys.argv[1:] or [os.path.basename(p).replace("_traj.dcd", "")
                            for p in sorted(glob.glob(os.path.join(MD, "*_traj.dcd")))
                            if "_rep" not in p]
    out = []
    for t in tags:
        try:
            r = run(t)
        except Exception as e:
            print(f"[fail] {t}: {e}")
            continue
        if not r:
            continue
        out.append(r)
        top5 = r["residues"][:5]
        print(f"[ok] {r['label']:<4} {t:<16} total {r['total_kJ_mol']:>9.1f} kJ/mol | "
              f"top: " + ", ".join(f"{z['resname']}{z['resid']}({z['energy_kJ_mol']:.0f})"
                                   for z in top5), flush=True)
    with open(os.path.join(OUT, "interface_energetics.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {os.path.join(OUT, 'interface_energetics.json')}")
