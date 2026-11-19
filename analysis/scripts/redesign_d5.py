"""
Hypothesis-driven redesign round for D5.

The escape map found D5's liability: it leans on gp120 Arg444 (+6.6 kcal/mol at a
position conserved in only 40% of Env sequences) and Lys305 (+8.3 kcal/mol, 57%).
Those are cheap escape routes. The virus cannot be mutated, so the fix must come
from the binder: reduce D5's dependence on those contacts while preserving overall
binding.

Procedure
  1. identify D5 residues contacting the liability positions
  2. build point mutants with PDBFixer (proper sidechain geometry)
  3. minimise, equilibrate briefly, and run short MD
  4. score with MM/GBSA and re-measure dependence on the liability contacts

A mutant is an improvement only if total binding holds AND dependence on the
variable positions falls. Losing binding to lose the liability is not progress,
and the script reports both so that trade-off is visible rather than hidden.
"""
import os, sys, json, warnings
import numpy as np
warnings.filterwarnings("ignore")

import mdtraj as md
from openmm import LangevinMiddleIntegrator, Platform, CustomExternalForce
from openmm.app import (ForceField, Simulation, PDBFile, HBonds, NoCutoff,
                        Modeller, DCDReporter)
from openmm.unit import (kelvin, picoseconds, femtoseconds, nanometer,
                         kilojoule_per_mole, amu)
from pdbfixer import PDBFixer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from alanine_scan import Decoupler, energy, BACKBONE

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SPLIT = os.path.join(PROJECT, "paper_analysis", "01_structures")
MD = os.path.join(PROJECT, "paper_analysis", "03_md")
OUT = os.path.join(PROJECT, "paper_analysis", "07_redesign")
os.makedirs(OUT, exist_ok=True)

BASE = "Latest_Best"
LIABILITY = [444, 305, 440, 308]      # gp120 positions, <60% conserved, load-bearing
NS = 0.4
FF = ("amber14-all.xml", "implicit/gbn2.xml")


def fix(path, mutations=None, chain="B"):
    fx = PDBFixer(filename=path)
    if mutations:
        fx.applyMutations(mutations, chain)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)
    return fx.topology, fx.positions


def make_sim(top, pos, restrain_distal=True):
    ff = ForceField(*FF)
    system = ff.createSystem(top, nonbondedMethod=NoCutoff, constraints=HBonds,
                             hydrogenMass=1.5 * amu,
                             soluteDielectric=1.0, solventDielectric=78.5)
    if restrain_distal:
        atoms = list(top.atoms())
        xyz = np.array(pos.value_in_unit(nanometer))
        lig = [i for i, a in enumerate(atoms)
               if a.residue.chain.id == "B" and a.element.symbol != "H"]
        L = xyz[lig]
        r = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        r.addGlobalParameter("k", 10.0 * kilojoule_per_mole / nanometer**2)
        for p in ("x0", "y0", "z0"):
            r.addPerParticleParameter(p)
        for i, a in enumerate(atoms):
            if a.residue.chain.id != "A" or a.name not in ("CA", "N", "C", "O"):
                continue
            if np.linalg.norm(L - xyz[i], axis=1).min() < 1.5:
                continue
            r.addParticle(i, xyz[i].tolist())
        system.addForce(r)
    integ = LangevinMiddleIntegrator(310.15 * kelvin, 1.0 / picoseconds, 4.0 * femtoseconds)
    try:
        sim = Simulation(top, system, integ, Platform.getPlatformByName("OpenCL"))
    except Exception:
        sim = Simulation(top, system, integ)
    sim.context.setPositions(pos)
    return sim


def hxb2_to_fixed():
    """Map HXB2 gp120 numbering to PDBFixer's output numbering.

    PDBFixer renumbers chain A sequentially from 1, so the HXB2 ids used by the
    escape map (Arg444, Lys305, ...) address different residues in the fixed
    topology. Residue order is preserved, so the mapping is positional.
    """
    orig = PDBFile(os.path.join(SPLIT, f"{BASE}_split.pdb")).topology
    ids = [int(r.id) for r in orig.residues() if r.chain.id == "A"]
    return {h: i + 1 for i, h in enumerate(ids)}


def score(tag, top, traj_xyz):
    """MM/GBSA plus dependence on the liability contacts."""
    atoms = list(top.atoms())
    rec = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "A"])
    lig = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "B"])
    mod_r = Modeller(top, [__import__("openmm").Vec3(0, 0, 0)] * top.getNumAtoms())
    mod_r.delete([r for r in top.residues() if r.chain.id != "A"])
    mod_l = Modeller(top, [__import__("openmm").Vec3(0, 0, 0)] * top.getNumAtoms())
    mod_l.delete([r for r in top.residues() if r.chain.id != "B"])

    p0 = [__import__("openmm").Vec3(*v) for v in traj_xyz[0]]
    sc = make_sim(top, p0, restrain_distal=False)
    sr = make_sim(mod_r.topology, [p0[i] for i in rec], restrain_distal=False)
    sl = make_sim(mod_l.topology, [p0[i] for i in lig], restrain_distal=False)

    X = traj_xyz
    e_rec = np.array([energy(sr, x[rec]) for x in X])
    wt = np.array([energy(sc, x) for x in X]) - e_rec - \
         np.array([energy(sl, x[lig]) for x in X])

    dec_c, dec_r = Decoupler(sc), Decoupler(sr)
    rec_pos = {g: k for k, g in enumerate(rec)}
    dep = {}
    remap = hxb2_to_fixed()
    for pos_id in LIABILITY:
        want = remap.get(pos_id, pos_id)
        res = next((a.residue for a in atoms
                    if a.residue.chain.id == "A" and int(a.residue.id) == want), None)
        if res is None:
            continue
        side = [a.index for a in res.atoms() if a.name not in BACKBONE]
        if not side:
            continue
        side_r = [rec_pos[i] for i in side if i in rec_pos]
        dec_c.zero(side); dec_r.zero(side_r)
        mut = np.array([energy(sc, x) for x in X]) - \
              np.array([energy(sr, x[rec]) for x in X]) - \
              np.array([energy(sl, x[lig]) for x in X])
        dec_c.restore(); dec_r.restore()
        dep[f"{res.name}{pos_id}"] = round(float((mut - wt).mean()) / 4.184, 2)
    sem = float(wt.std(ddof=1) / np.sqrt(len(wt))) if len(wt) > 1 else 0.0
    return float(wt.mean()), sem, dep


def run_variant(name, mutations):
    src = os.path.join(SPLIT, f"{BASE}_split.pdb")
    top, pos = fix(src, mutations)
    sim = make_sim(top, pos)
    sim.minimizeEnergy(maxIterations=1500)
    sim.context.setVelocitiesToTemperature(310.15 * kelvin)
    sim.step(int(50 * 1000 / 4.0))                       # 50 ps equilibration
    dcd = os.path.join(OUT, f"{name}_traj.dcd")
    sim.reporters.append(DCDReporter(dcd, 2500, enforcePeriodicBox=False))
    sim.step(int(NS * 1e6 / 4.0))
    # On Windows the DCDReporter keeps its file handle open; mdtraj then cannot
    # open the trajectory ("premature end of file") and the whole run is wasted.
    for rep in sim.reporters:
        try:
            rep._out.close()
        except Exception:
            pass
    sim.reporters.clear()

    refpdb = os.path.join(OUT, f"{name}_top.pdb")
    with open(refpdb, "w") as fh:
        PDBFile.writeFile(top, pos, fh)
    tr = md.load(dcd, top=refpdb)
    keep = tr.xyz[len(tr) // 2:]
    dg, sem, dep = score(name, top, keep)
    return {"variant": name, "mutations": mutations,
            "dG_kJ_mol": round(dg, 1), "dG_kcal_mol": round(dg / 4.184, 1),
            "sem_kcal_mol": round(sem / 4.184, 1),
            "liability_dependence_kcal": dep,
            "liability_total_kcal": round(sum(dep.values()), 2),
            "n_frames": len(keep)}


VARIANTS = [
    ("WT", []),
    ("R40A", ["ARG-40-ALA"]),
    ("K39A", ["LYS-39-ALA"]),
    ("E62Q", ["GLU-62-GLN"]),
    ("R40E", ["ARG-40-GLU"]),
    ("K66A", ["LYS-66-ALA"]),
    ("R40A_K39A", ["ARG-40-ALA", "LYS-39-ALA"]),
]

if __name__ == "__main__":
    out = []
    for name, muts in VARIANTS:
        try:
            r = run_variant(name, muts)
        except Exception as e:
            print(f"[fail] {name}: {e}", flush=True)
            continue
        out.append(r)
        print(f"[ok] {name:<12} dG {r['dG_kcal_mol']:>8.1f} kcal/mol | "
              f"liability {r['liability_total_kcal']:>6.2f} | {r['liability_dependence_kcal']}",
              flush=True)
        with open(os.path.join(OUT, "redesign_d5.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("REDESIGN DONE", flush=True)
