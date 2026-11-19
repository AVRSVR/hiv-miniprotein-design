"""
MM/GBSA binding free energies from the MD trajectories.

Single-trajectory MM/GBSA: for every saved frame, the potential energy of the
complex, the receptor alone and the binder alone are evaluated with the same
force field and implicit-solvent model, and

    dG_bind = <E_complex> - <E_receptor> - <E_binder>

Because receptor and binder coordinates are taken from the complex trajectory,
intramolecular strain largely cancels and the result is dominated by interaction
and desolvation terms.

Interpretation limits - these matter and are stated in the manuscript:
  * no configurational entropy term (no normal-mode or quasi-harmonic correction),
    so these are enthalpy-dominated effective energies, not absolute affinities;
  * single-trajectory, so conformational adaptation on binding is not captured;
  * implicit solvent omits explicit interfacial water and ion effects.
MM/GBSA of this form is useful for *ranking* related complexes, not for predicting
Kd. Values are reported with the standard error over frames, and the per-frame
spread is large by construction - report it.

Usage: python mmgbsa.py [tag ...]
"""
import os, sys, json, glob, warnings
import numpy as np
warnings.filterwarnings("ignore")

import mdtraj as md
from openmm import Platform, LangevinMiddleIntegrator
from openmm.app import ForceField, Modeller, HBonds, NoCutoff, PDBFile, Simulation
from openmm.unit import kelvin, picoseconds, femtoseconds, nanometer, kilojoule_per_mole, amu
from pdbfixer import PDBFixer

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SPLIT = os.path.join(PROJECT, "paper_analysis", "01_structures")
MD = os.path.join(PROJECT, "paper_analysis", "03_md")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")

KJ = kilojoule_per_mole
FF = ("amber14-all.xml", "implicit/gbn2.xml")


def build_topology(split_pdb):
    """Reproduce exactly the topology run_md.py simulated."""
    fx = PDBFixer(filename=split_pdb)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)
    return fx.topology, fx.positions


def make_context(topology, positions):
    ff = ForceField(*FF)
    system = ff.createSystem(topology, nonbondedMethod=NoCutoff, constraints=HBonds,
                             hydrogenMass=1.5 * amu,
                             soluteDielectric=1.0, solventDielectric=78.5)
    integ = LangevinMiddleIntegrator(310.15 * kelvin, 1.0 / picoseconds, 4.0 * femtoseconds)
    try:
        sim = Simulation(topology, system, integ, Platform.getPlatformByName("OpenCL"))
    except Exception:
        sim = Simulation(topology, system, integ)
    sim.context.setPositions(positions)
    return sim


def subset(topology, positions, keep_chain):
    mod = Modeller(topology, positions)
    mod.delete([r for r in topology.residues() if r.chain.id != keep_chain])
    return mod.topology, mod.positions


def energy(sim, xyz_nm):
    sim.context.setPositions(xyz_nm)
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)


def run(tag):
    split_pdb = os.path.join(SPLIT, f"{tag}_split.pdb")
    dcd = os.path.join(MD, f"{tag}_traj.dcd")
    if not (os.path.exists(split_pdb) and os.path.exists(dcd)):
        print(f"[skip] {tag}: missing inputs")
        return None

    top, pos = build_topology(split_pdb)
    ref_pdb = os.path.join(MD, f"{tag}_mmgbsa_top.pdb")
    with open(ref_pdb, "w") as fh:
        PDBFile.writeFile(top, pos, fh)

    traj = md.load(dcd, top=ref_pdb)
    if traj.n_frames < 5:
        print(f"[skip] {tag}: only {traj.n_frames} frames")
        return None

    atoms = list(top.atoms())
    rec_idx = [i for i, a in enumerate(atoms) if a.residue.chain.id == "A"]
    lig_idx = [i for i, a in enumerate(atoms) if a.residue.chain.id == "B"]
    if not rec_idx or not lig_idx:
        print(f"[skip] {tag}: chains not found")
        return None

    top_r, pos_r = subset(top, pos, "A")
    top_l, pos_l = subset(top, pos, "B")

    sim_c = make_context(top, pos)
    sim_r = make_context(top_r, pos_r)
    sim_l = make_context(top_l, pos_l)

    # discard the first 25% as equilibration within the production segment
    start = max(1, traj.n_frames // 4)
    dgs = []
    for f in range(start, traj.n_frames):
        x = traj.xyz[f]
        try:
            ec = energy(sim_c, x)
            er = energy(sim_r, x[rec_idx])
            el = energy(sim_l, x[lig_idx])
        except Exception:
            continue
        dgs.append(ec - er - el)

    if len(dgs) < 3:
        print(f"[skip] {tag}: too few usable frames")
        return None

    d = np.array(dgs)
    res = {"design": tag, "n_frames": len(d),
           "dG_bind_kJ_mol": round(float(d.mean()), 1),
           "sem_kJ_mol": round(float(d.std(ddof=1) / np.sqrt(len(d))), 1),
           "sd_kJ_mol": round(float(d.std(ddof=1)), 1),
           "dG_bind_kcal_mol": round(float(d.mean()) / 4.184, 1),
           "sem_kcal_mol": round(float(d.std(ddof=1) / np.sqrt(len(d))) / 4.184, 1)}
    print(f"[ok] {tag:<16} dG = {res['dG_bind_kcal_mol']:8.1f} +/- "
          f"{res['sem_kcal_mol']:.1f} kcal/mol  ({len(d)} frames)")
    return res


if __name__ == "__main__":
    tags = sys.argv[1:]
    if not tags:
        tags = [os.path.basename(p).replace("_traj.dcd", "")
                for p in sorted(glob.glob(os.path.join(MD, "*_traj.dcd")))]
    rows = [r for r in (run(t) for t in tags) if r]
    if rows:
        with open(os.path.join(OUT, "mmgbsa.json"), "w") as f:
            json.dump(rows, f, indent=2)
        print(f"\nwrote {os.path.join(OUT, 'mmgbsa.json')}")
