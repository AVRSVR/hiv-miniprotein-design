"""
Computational alanine scanning of the binder interface, by sidechain decoupling.

Per-residue decomposition of a raw interaction energy is dominated by charged
residues: the Coulomb gain of a salt bridge is not offset by the desolvation
penalty of burying the charge. Alanine scanning avoids this, because the GB term
supplies the missing desolvation, and

    ddG = dG(mutant) - dG(wild type)

is positive for residues whose removal costs binding - the hotspots.

The mutant is produced by *decoupling* rather than deletion: for every sidechain
atom beyond CB, the partial charge and Lennard-Jones well depth are set to zero in
both the NonbondedForce and the generalised-Born force. The atoms remain present
but cease to interact, so atom count and ordering are unchanged and the
wild-type trajectory frames map onto the mutant without rebuilding topology or
re-adding hydrogens. Parameters are pushed into the existing contexts with
updateParametersInContext, so the whole scan reuses three contexts.

Two consequences to read the numbers with: the decoupled sidechain still occupies
volume (it simply does not interact), and no repacking or relaxation follows the
mutation. ddG is therefore an estimate of the direct nonbonded contribution of the
sidechain, not the free energy of a real point mutant.

Usage: python alanine_scan.py <tag> [n_frames]
"""
import os, sys, json, warnings
import numpy as np
warnings.filterwarnings("ignore")

import mdtraj as md
from openmm import Platform, LangevinMiddleIntegrator, NonbondedForce, CustomGBForce, GBSAOBCForce
from openmm.app import ForceField, Modeller, Simulation, PDBFile, HBonds, NoCutoff
from openmm.unit import kelvin, picoseconds, femtoseconds, kilojoule_per_mole, amu
from pdbfixer import PDBFixer

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SPLIT = os.path.join(PROJECT, "paper_analysis", "01_structures")
MD = os.path.join(PROJECT, "paper_analysis", "03_md")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")

BACKBONE = {"N", "CA", "C", "O", "CB", "OXT", "H", "HA", "HA2", "HA3",
            "HB", "HB1", "HB2", "HB3"}
SKIP_RES = {"GLY", "ALA", "PRO"}
IFACE_CUT = 0.45
FF = ("amber14-all.xml", "implicit/gbn2.xml")
LABEL = {"S310": "D1", "Chai1_Best": "D2", "Cyclic_Ring": "D3",
         "model_4_chai": "D4", "Latest_Best": "D5"}


def build_wt(split_pdb):
    fx = PDBFixer(filename=split_pdb)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)
    return fx.topology, fx.positions


def make_sim(topology, positions):
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


class Decoupler:
    """Zero / restore nonbonded + GB parameters for a set of atom indices."""

    def __init__(self, sim):
        self.sim = sim
        self.nb = None
        self.gb = None
        for f in sim.system.getForces():
            if isinstance(f, NonbondedForce):
                self.nb = f
            elif isinstance(f, (CustomGBForce, GBSAOBCForce)):
                self.gb = f
        self.saved = {}

    def zero(self, idxs):
        self.saved = {"nb": {}, "gb": {}}
        for i in idxs:
            q, s, e = self.nb.getParticleParameters(i)
            self.saved["nb"][i] = (q, s, e)
            self.nb.setParticleParameters(i, 0.0, s, 0.0)
        self.nb.updateParametersInContext(self.sim.context)
        if self.gb is not None:
            for i in idxs:
                p = list(self.gb.getParticleParameters(i))
                self.saved["gb"][i] = list(p)
                p[0] = 0.0                      # charge
                self.gb.setParticleParameters(i, p)
            self.gb.updateParametersInContext(self.sim.context)

    def restore(self):
        for i, (q, s, e) in self.saved.get("nb", {}).items():
            self.nb.setParticleParameters(i, q, s, e)
        self.nb.updateParametersInContext(self.sim.context)
        if self.gb is not None:
            for i, p in self.saved.get("gb", {}).items():
                self.gb.setParticleParameters(i, p)
            self.gb.updateParametersInContext(self.sim.context)
        self.saved = {}


def energy(sim, x):
    sim.context.setPositions(x)
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        kilojoule_per_mole)


def main(tag, n_frames=20):
    split_pdb = os.path.join(SPLIT, f"{tag}_split.pdb")
    dcd = os.path.join(MD, f"{tag}_traj.dcd")
    top, pos = build_wt(split_pdb)
    ref = os.path.join(MD, f"{tag}_mmgbsa_top.pdb")
    if not os.path.exists(ref):
        with open(ref, "w") as fh:
            PDBFile.writeFile(top, pos, fh)
    traj = md.load(dcd, top=ref)

    atoms = list(top.atoms())
    rec_i = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "A"])
    lig_i = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "B"])

    idx = np.linspace(traj.n_frames // 2, traj.n_frames - 1, n_frames).astype(int)
    X = traj.xyz[idx]

    # ligand-only topology, positional mapping preserved by Modeller.delete
    mod = Modeller(top, pos)
    mod.delete([r for r in top.residues() if r.chain.id != "B"])
    top_l, pos_l = mod.topology, mod.positions
    lig_map = lig_i                                  # complex index -> ligand order
    lig_pos_in_lig = {g: k for k, g in enumerate(lig_i)}

    sim_c = make_sim(top, pos)
    sim_l = make_sim(top_l, pos_l)
    mod_r = Modeller(top, pos)
    mod_r.delete([r for r in top.residues() if r.chain.id != "A"])
    sim_r = make_sim(mod_r.topology, mod_r.positions)

    # receptor is never mutated, so its energies are computed once
    e_rec = np.array([energy(sim_r, x[rec_i]) for x in X])

    dec_c = Decoupler(sim_c)
    dec_l = Decoupler(sim_l)

    wt = np.array([energy(sim_c, x) for x in X]) - e_rec - \
         np.array([energy(sim_l, x[lig_map]) for x in X])
    print(f"{tag}: wild-type dG = {wt.mean():.1f} kJ/mol over {len(X)} frames", flush=True)

    # interface residues on the binder
    x0 = X[len(X) // 2]
    heavy_r = [i for i in rec_i if atoms[i].element.symbol != "H"]
    iface = set()
    for i in lig_i:
        if atoms[i].element.symbol == "H":
            continue
        if np.linalg.norm(x0[heavy_r] - x0[i], axis=1).min() < IFACE_CUT:
            iface.add(atoms[i].residue.index)
    targets = [r for r in sorted(iface)
               if next(a for a in atoms if a.residue.index == r).residue.name not in SKIP_RES]
    print(f"  {len(iface)} interface residues, {len(targets)} scannable", flush=True)

    rows = []
    for n, rid in enumerate(targets):
        res = next(a.residue for a in atoms if a.residue.index == rid)
        side_c = [a.index for a in res.atoms() if a.name not in BACKBONE]
        if not side_c:
            continue
        side_l = [lig_pos_in_lig[i] for i in side_c if i in lig_pos_in_lig]

        dec_c.zero(side_c)
        dec_l.zero(side_l)
        mut = np.array([energy(sim_c, x) for x in X]) - e_rec - \
              np.array([energy(sim_l, x[lig_map]) for x in X])
        dec_c.restore()
        dec_l.restore()

        d = mut - wt
        ddg = float(d.mean())
        rows.append({"resname": res.name, "resid": int(res.id),
                     "ddG_kJ_mol": round(ddg, 2),
                     "ddG_kcal_mol": round(ddg / 4.184, 2),
                     "sem_kcal_mol": round(float(np.std(d, ddof=1) /
                                                 np.sqrt(len(d)) / 4.184), 2),
                     "hotspot": bool(ddg >= 4.184)})
        print(f"  [{n+1}/{len(targets)}] {res.name}{res.id}: "
              f"{ddg/4.184:+.2f} kcal/mol", flush=True)

    rows.sort(key=lambda z: -z["ddG_kJ_mol"])
    out = {"design": tag, "label": LABEL.get(tag, tag), "n_frames": len(X),
           "wt_dG_kJ_mol": round(float(wt.mean()), 1),
           "n_hotspots": sum(1 for z in rows if z["hotspot"]), "residues": rows}
    path = os.path.join(OUT, f"alascan_{tag}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    hs = [z for z in rows if z["hotspot"]]
    print(f"\n{tag}: {len(hs)} hotspots >=1 kcal/mol: " +
          ", ".join(f"{z['resname']}{z['resid']}({z['ddG_kcal_mol']:+.1f})" for z in hs[:12]))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 20)
