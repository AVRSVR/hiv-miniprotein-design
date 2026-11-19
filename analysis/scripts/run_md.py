"""
Interface-stability MD for designed miniprotein / gp120 complexes.

System setup
------------
GBn2 implicit solvent (amber14 + implicit/gbn2.xml), 310.15 K Langevin, 2 fs with
HBonds constrained. gp120 backbone atoms further than RESTRAIN_CUT from the binder
are positionally restrained: gp120 is natively trimer-embedded and glycosylated, so
restraining the distal scaffold prevents implicit-solvent unfolding artifacts while
leaving the interface free to relax.

Per-frame metrics
-----------------
  iRMSD        binder CA deviation after superposing on receptor interface CA
  binder_rmsd  binder internal fold deviation (self-aligned)
  contacts     receptor-binder heavy-atom pairs < 4.5 A
  binder_rg    binder radius of gyration

Usage: python run_md.py <split_pdb> <tag> <ns>
"""
import os, sys, json, time
import numpy as np
from openmm import (LangevinMiddleIntegrator, Platform, CustomExternalForce, unit)
from openmm.app import (PDBFile, ForceField, Simulation, DCDReporter, HBonds, NoCutoff)
from openmm.unit import kelvin, picoseconds, femtoseconds, nanometer, kilojoule_per_mole, amu
from pdbfixer import PDBFixer

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
OUT = os.path.join(PROJECT, "paper_analysis", "03_md")
os.makedirs(OUT, exist_ok=True)

TEMP = 310.15
DT = 4.0                    # fs, enabled by hydrogen mass repartitioning
REPORT_PS = 10.0            # frame every 10 ps
EQUIL_PS = 100.0            # ps equilibration
RESTRAIN_CUT = 1.5          # nm from binder; beyond this gp120 backbone is restrained
RESTRAIN_K = 10.0           # kJ/mol/nm^2
CONTACT_CUT = 0.45          # nm


def kabsch(P, Q):
    """Return P superposed onto Q (both (N,3))."""
    pc, qc = P.mean(0), Q.mean(0)
    p, q = P - pc, Q - qc
    V, S, Wt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    return (p @ R) + qc


def rmsd_after_fit(mob, ref, fit_mob, fit_ref):
    """RMSD of mob vs ref using the transform that fits fit_mob onto fit_ref."""
    pc, qc = fit_mob.mean(0), fit_ref.mean(0)
    p, q = fit_mob - pc, fit_ref - qc
    V, S, Wt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    moved = ((mob - pc) @ R) + qc
    return float(np.sqrt(((moved - ref) ** 2).sum(1).mean()))


def main(split_pdb, tag, ns):
    t0 = time.time()
    print(f"=== MD: {tag} ({ns} ns) ===", flush=True)

    fixer = PDBFixer(filename=split_pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # don't build long absent loops de novo
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    top, pos = fixer.topology, fixer.positions
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(top, nonbondedMethod=NoCutoff,
                             constraints=HBonds, hydrogenMass=1.5*amu,
                             soluteDielectric=1.0, solventDielectric=78.5)

    atoms = list(top.atoms())
    chain_ids = [a.residue.chain.id for a in atoms]
    xyz = np.array(pos.value_in_unit(nanometer))

    rec_idx = [i for i, a in enumerate(atoms) if chain_ids[i] == "A"]
    lig_idx = [i for i, a in enumerate(atoms) if chain_ids[i] == "B"]
    if not rec_idx or not lig_idx:
        raise SystemExit("need chains A (receptor) and B (binder)")

    rec_heavy = [i for i in rec_idx if atoms[i].element is not None and atoms[i].element.symbol != "H"]
    lig_heavy = [i for i in lig_idx if atoms[i].element is not None and atoms[i].element.symbol != "H"]
    rec_ca = [i for i in rec_idx if atoms[i].name == "CA"]
    lig_ca = [i for i in lig_idx if atoms[i].name == "CA"]

    # distance of each receptor atom to nearest binder atom (initial pose)
    lig_xyz0 = xyz[lig_heavy]
    d_rec = np.linalg.norm(xyz[rec_heavy][:, None, :] - lig_xyz0[None, :, :], axis=-1).min(1)
    near_rec_heavy = set(np.array(rec_heavy)[d_rec < RESTRAIN_CUT].tolist())

    # interface receptor CA (for iRMSD superposition): CA of residues with a near atom
    near_res = {atoms[i].residue.index for i in near_rec_heavy}
    iface_rec_ca = [i for i in rec_ca if atoms[i].residue.index in near_res]
    if len(iface_rec_ca) < 4:
        iface_rec_ca = rec_ca
    print(f"receptor atoms {len(rec_idx)}  binder atoms {len(lig_idx)}  "
          f"iface rec CA {len(iface_rec_ca)}", flush=True)

    # restrain distal receptor backbone
    restraint = CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    restraint.addGlobalParameter("k", RESTRAIN_K * kilojoule_per_mole / nanometer**2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    n_restr = 0
    for i in rec_idx:
        if atoms[i].name not in ("CA", "N", "C", "O"):
            continue
        if i in near_rec_heavy:
            continue
        restraint.addParticle(i, xyz[i].tolist())
        n_restr += 1
    system.addForce(restraint)
    print(f"restrained {n_restr} distal receptor backbone atoms", flush=True)

    integrator = LangevinMiddleIntegrator(TEMP * kelvin, 1.0 / picoseconds, DT * femtoseconds)
    try:
        plat = Platform.getPlatformByName("OpenCL")
        sim = Simulation(top, system, integrator, plat)
    except Exception:
        sim = Simulation(top, system, integrator)
    print("platform:", sim.context.getPlatform().getName(), flush=True)

    sim.context.setPositions(pos)
    print("minimizing...", flush=True)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(TEMP * kelvin)

    equil_steps = int(EQUIL_PS * 1000 / DT)
    print(f"equilibrating {EQUIL_PS} ps...", flush=True)
    sim.step(equil_steps)

    ref = np.array(sim.context.getState(getPositions=True)
                   .getPositions().value_in_unit(nanometer))
    ref_lig_ca, ref_iface_ca = ref[lig_ca], ref[iface_rec_ca]

    dcd = os.path.join(OUT, f"{tag}_traj.dcd")
    # The system is non-periodic (implicit solvent, NoCutoff) but still carries
    # OpenMM's default 2 nm box. Left to its default, DCDReporter wraps every
    # molecule into that box and the ~10 nm complex is torn apart in the saved
    # file - the simulation is unaffected, but the trajectory becomes unusable.
    sim.reporters.append(DCDReporter(dcd, int(REPORT_PS * 1000 / DT),
                                     enforcePeriodicBox=False))

    total_steps = int(ns * 1e6 / DT)
    per_frame = int(REPORT_PS * 1000 / DT)
    n_frames = total_steps // per_frame

    rows = []
    for f in range(n_frames):
        sim.step(per_frame)
        st = sim.context.getState(getPositions=True, getEnergy=True)
        c = np.array(st.getPositions().value_in_unit(nanometer))

        irmsd = rmsd_after_fit(c[lig_ca], ref_lig_ca, c[iface_rec_ca], ref_iface_ca) * 10.0
        b = c[lig_ca]
        brmsd = float(np.sqrt(((kabsch(b, ref_lig_ca) - ref_lig_ca) ** 2).sum(1).mean())) * 10.0
        dm = np.linalg.norm(c[rec_heavy][:, None, :] - c[lig_heavy][None, :, :], axis=-1)
        contacts = int((dm < CONTACT_CUT).sum())
        mind = float(dm.min()) * 10.0
        bc = b - b.mean(0)
        rg = float(np.sqrt((bc ** 2).sum(1).mean())) * 10.0

        rows.append({"time_ps": (f + 1) * REPORT_PS, "iRMSD_A": irmsd,
                     "binder_RMSD_A": brmsd, "contacts": contacts,
                     "min_dist_A": mind, "binder_Rg_A": rg,
                     "PE_kJ": st.getPotentialEnergy().value_in_unit(kilojoule_per_mole)})

        if (f + 1) % 25 == 0 or f == n_frames - 1:
            el = time.time() - t0
            print(f"  {(f+1)*REPORT_PS/1000:.2f} ns | iRMSD {irmsd:.2f} A | "
                  f"contacts {contacts} | mind {mind:.2f} A | {el/60:.1f} min", flush=True)
            # checkpoint: a killed run should still leave usable metrics
            with open(os.path.join(OUT, f"{tag}_md_metrics.json"), "w") as fh:
                json.dump({"design": tag, "ns": (f + 1) * REPORT_PS / 1000.0,
                           "ns_target": ns, "complete": f == n_frames - 1,
                           "temp_K": TEMP, "solvent": "GBn2 implicit",
                           "frames": rows}, fh, indent=2)

    with open(os.path.join(OUT, f"{tag}_md_metrics.json"), "w") as fh:
        json.dump({"design": tag, "ns": ns, "temp_K": TEMP,
                   "solvent": "GBn2 implicit", "frames": rows}, fh, indent=2)

    # final frame for figures
    with open(os.path.join(OUT, f"{tag}_final.pdb"), "w") as fh:
        PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(), fh)

    print(f"done {tag} in {(time.time()-t0)/60:.1f} min -> {dcd}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
