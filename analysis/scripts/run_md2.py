"""
Interface-stability MD for designed miniprotein / gp120 complexes.

Setup
-----
The gp120 scaffold is truncated to the interface neighbourhood (residues with any
atom within TRUNC_CUT of the binder in the starting pose). gp120 is natively
trimer-embedded and glycosylated, so the distal scaffold is not free to move in
reality; keeping it would only add cost. Residues in the outer shell of the
retained region are positionally restrained to hold the epitope geometry, leaving
the interface itself free to relax.

GBn2 implicit solvent, 310.15 K Langevin, hydrogen mass repartitioning for a 4 fs
timestep with HBonds constrained.

Per-frame metrics
-----------------
  iRMSD        binder CA deviation after superposing on receptor interface CA
  binder_RMSD  binder internal fold deviation
  contacts     receptor-binder heavy-atom pairs < 4.5 A
  min_dist     closest receptor-binder heavy-atom approach
  binder_Rg    binder radius of gyration

Usage: python run_md2.py <split_pdb> <tag> <ns>
"""
import os, sys, json, time
import numpy as np
from openmm import LangevinMiddleIntegrator, Platform, CustomExternalForce
from openmm.app import PDBFile, ForceField, Simulation, DCDReporter, HBonds, NoCutoff, Modeller
from openmm.unit import (kelvin, picoseconds, femtoseconds, nanometer,
                         kilojoule_per_mole, amu, dalton)
from pdbfixer import PDBFixer

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
OUT = os.path.join(PROJECT, "paper_analysis", "03_md")
os.makedirs(OUT, exist_ok=True)

TEMP = 310.15
DT = float(os.environ.get("MD_DT_FS", "4.0"))   # fs, enabled by hydrogen mass repartitioning
HMASS = float(os.environ.get("MD_HMASS", "1.5"))  # amu
REPORT_PS = 10.0
EQUIL_PS = 100.0
TRUNC_CUT = 1.2          # nm: receptor residues kept
SHELL_RESTRAIN = 0.9     # nm: kept residues beyond this are restrained
MIN_FRAG = 3             # residues: shorter receptor fragments are discarded
# Defaults reproduce the published panel runs exactly. The tandem construct needs
# a gentler start because its linker is built by hand, so those extras are opt-in
# through the environment rather than changed for everything.
MIN_ITERS = int(os.environ.get("MD_MIN_ITERS", "1500"))
RAMP = float(os.environ.get("MD_RAMP_PS", "0"))  # ps per heating stage; 0 disables
RESTRAIN_K = 50.0
CONTACT_CUT = 0.45


def _write_fragmented(topology, positions, out_pdb, min_frag=MIN_FRAG):
    """
    Write a PDB where each contiguous receptor fragment becomes its own chain, so
    PDBFixer will cap every cut end as a real terminus. The binder keeps chain 'B'.

    Fragments shorter than min_frag residues are discarded. A singleton gets both
    N-terminal hydrogens and an OXT from the capping step, i.e. a free zwitterionic
    amino acid, for which amber14 has no template and createSystem fails. Such
    fragments only appear when the binder is large enough for the truncation
    envelope to catch scattered receptor residues; they sit in the restrained outer
    shell, far outside contact range, so dropping them does not touch the interface.
    Dropped residues are listed so the omission is on the record.
    """
    import string
    pos = np.array(positions.value_in_unit(nanometer))
    pool = [c for c in (string.ascii_uppercase + string.digits) if c != "B"]

    rec_res = [r for r in topology.residues() if r.chain.id != "B"]
    lig_res = [r for r in topology.residues() if r.chain.id == "B"]

    # split receptor into fragments on backbone discontinuity (C_i -> N_i+1)
    frags, cur = [], []
    for r in rec_res:
        if not cur:
            cur = [r]
            continue
        prev = cur[-1]
        c_at = next((a for a in prev.atoms() if a.name == "C"), None)
        n_at = next((a for a in r.atoms() if a.name == "N"), None)
        joined = (c_at is not None and n_at is not None and
                  np.linalg.norm(pos[c_at.index] - pos[n_at.index]) < 0.2)
        if joined:
            cur.append(r)
        else:
            frags.append(cur); cur = [r]
    if cur:
        frags.append(cur)

    dropped = [r for f in frags if len(f) < min_frag for r in f]
    if dropped:
        lig_xyz = np.array([pos[a.index] for r in lig_res for a in r.atoms()
                            if a.element is not None and a.element.symbol != "H"])
        tags = []
        for r in dropped:
            rx = np.array([pos[a.index] for a in r.atoms()
                           if a.element is not None and a.element.symbol != "H"])
            d = np.linalg.norm(rx[:, None, :] - lig_xyz[None, :, :], axis=-1).min()
            tags.append("%s%s (%.2f nm from binder)" % (r.name, r.id, d))
        print("  dropping %d residue(s) in fragments shorter than %d: %s"
              % (len(dropped), min_frag, ", ".join(tags)))
    frags = [f for f in frags if len(f) >= min_frag]

    lines, serial = [], 1
    for fi, frag in enumerate(frags):
        cid = pool[fi % len(pool)]
        for r in frag:
            for a in r.atoms():
                x, y, z = pos[a.index] * 10.0
                el = a.element.symbol if a.element is not None else " "
                lines.append(
                    f"ATOM  {serial:5d} {a.name:^4s} {r.name:>3s} {cid}{r.id:>4s}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
                serial += 1
        lines.append("TER\n")
    for r in lig_res:
        for a in r.atoms():
            x, y, z = pos[a.index] * 10.0
            el = a.element.symbol if a.element is not None else " "
            lines.append(
                f"ATOM  {serial:5d} {a.name:^4s} {r.name:>3s} B{r.id:>4s}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
            serial += 1
    lines.append("TER\nEND\n")
    with open(out_pdb, "w") as f:
        f.writelines(lines)
    return len(frags)


def kabsch_rmsd(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    p, q = P - pc, Q - qc
    V, S, Wt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return float(np.sqrt((((p @ R) - q) ** 2).sum(1).mean()))


def rmsd_after_fit(mob, ref, fit_mob, fit_ref):
    """RMSD of mob vs ref, using the rotation that superposes fit_mob onto fit_ref."""
    pc, qc = fit_mob.mean(0), fit_ref.mean(0)
    p, q = fit_mob - pc, fit_ref - qc
    V, S, Wt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    moved = ((mob - pc) @ R) + qc
    return float(np.sqrt(((moved - ref) ** 2).sum(1).mean()))


def main(split_pdb, tag, ns):
    t0 = time.time()
    log = open(os.path.join(OUT, f"{tag}_md.log"), "w", buffering=1)

    def say(m):
        print(m, flush=True)
        log.write(m + "\n")

    say(f"=== MD {tag} | target {ns} ns ===")

    # ---- truncate receptor around the binder -------------------------------
    fx = PDBFixer(filename=split_pdb)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)

    top0, pos0 = fx.topology, fx.positions
    xyz0 = np.array(pos0.value_in_unit(nanometer))
    atoms0 = list(top0.atoms())

    lig_heavy = [i for i, a in enumerate(atoms0)
                 if a.residue.chain.id == "B" and a.element is not None and a.element.symbol != "H"]
    if not lig_heavy:
        raise SystemExit("no binder chain B")
    L = xyz0[lig_heavy]

    keep_res, shell_res = set(), set()
    for res in top0.residues():
        if res.chain.id == "B":
            keep_res.add(res.index)
            continue
        ridx = [a.index for a in res.atoms()
                if a.element is not None and a.element.symbol != "H"]
        if not ridx:
            continue
        d = np.linalg.norm(xyz0[ridx][:, None, :] - L[None, :, :], axis=-1).min()
        if d < TRUNC_CUT:
            keep_res.add(res.index)
            if d > SHELL_RESTRAIN:
                shell_res.add(res.index)

    mod = Modeller(top0, pos0)
    mod.delete([r for r in top0.residues() if r.index not in keep_res])
    say(f"truncated receptor: kept {len(keep_res)} residues, "
        f"{mod.topology.getNumAtoms()} atoms (from {len(atoms0)})")

    # Truncation cuts the backbone, leaving fragment ends that are neither real
    # termini nor bonded to anything. Write each contiguous fragment as its own
    # chain, then re-fix so PDBFixer caps every fragment with proper termini.
    frag_pdb = os.path.join(OUT, f"{tag}_trunc.pdb")
    _write_fragmented(mod.topology, mod.positions, frag_pdb)

    fx2 = PDBFixer(filename=frag_pdb)
    fx2.findMissingResidues(); fx2.missingResidues = {}
    fx2.findMissingAtoms(); fx2.addMissingAtoms()
    fx2.addMissingHydrogens(7.4)
    top, pos = fx2.topology, fx2.positions
    atoms = list(top.atoms())
    xyz = np.array(pos.value_in_unit(nanometer))
    say(f"after capping: {len(atoms)} atoms, "
        f"{len(list(top.chains()))} chains, {len(list(top.residues()))} residues")

    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(top, nonbondedMethod=NoCutoff, constraints=HBonds,
                             hydrogenMass=HMASS * amu,
                             soluteDielectric=1.0, solventDielectric=78.5)

    # binder keeps chain B; receptor fragments occupy every other chain
    rec = [i for i, a in enumerate(atoms) if a.residue.chain.id != "B"]
    lig = [i for i, a in enumerate(atoms) if a.residue.chain.id == "B"]
    rec_heavy = [i for i in rec if atoms[i].element.symbol != "H"]
    lig_heavy2 = [i for i in lig if atoms[i].element.symbol != "H"]
    lig_ca = [i for i in lig if atoms[i].name == "CA"]

    # recompute distances in the rebuilt topology (indices changed on re-fix)
    Lh = xyz[lig_heavy2]
    d_rec = np.linalg.norm(xyz[rec_heavy][:, None, :] - Lh[None, :, :], axis=-1).min(1)
    near_res = {atoms[i].residue.index for i in np.array(rec_heavy)[d_rec < 0.5]}
    shell_now = {atoms[i].residue.index
                 for i in np.array(rec_heavy)[d_rec > SHELL_RESTRAIN]}
    iface_rec_ca = [i for i in rec if atoms[i].name == "CA"
                    and atoms[i].residue.index in near_res]
    if len(iface_rec_ca) < 4:
        iface_rec_ca = [i for i in rec if atoms[i].name == "CA"]

    restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", RESTRAIN_K * kilojoule_per_mole / nanometer**2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    n_r = 0
    for i in rec:
        if atoms[i].name not in ("CA", "N", "C", "O"):
            continue
        if atoms[i].residue.index in shell_now:
            restraint.addParticle(i, xyz[i].tolist())
            n_r += 1
    system.addForce(restraint)
    say(f"restrained {n_r} shell backbone atoms | iface rec CA {len(iface_rec_ca)} | "
        f"binder CA {len(lig_ca)}")

    integ = LangevinMiddleIntegrator(TEMP * kelvin, 1.0 / picoseconds, DT * femtoseconds)
    try:
        sim = Simulation(top, system, integ, Platform.getPlatformByName("OpenCL"))
    except Exception:
        sim = Simulation(top, system, integ)
    say(f"platform: {sim.context.getPlatform().getName()}")

    sim.context.setPositions(pos)
    sim.minimizeEnergy(maxIterations=MIN_ITERS)
    if RAMP:
        # A hand-built linker can arrive with strained backbone geometry, and
        # assigning 310 K velocities to it at a 4 fs timestep releases that strain
        # fast enough to integrate to NaN. Heating in stages lets it bleed off.
        for t_k in (50.0, 150.0, 250.0, TEMP):
            integ.setTemperature(t_k * kelvin)
            sim.context.setVelocitiesToTemperature(t_k * kelvin)
            sim.step(int(RAMP * 1000 / DT))
            sim.minimizeEnergy(maxIterations=200)
        say(f"ramped to {TEMP} K in 4 stages of {RAMP} ps")
    sim.context.setVelocitiesToTemperature(TEMP * kelvin)
    sim.step(int(EQUIL_PS * 1000 / DT))
    say(f"equilibrated {EQUIL_PS} ps ({(time.time()-t0)/60:.1f} min elapsed)")

    ref = np.array(sim.context.getState(getPositions=True)
                   .getPositions().value_in_unit(nanometer))
    ref_lig_ca, ref_if_ca = ref[lig_ca], ref[iface_rec_ca]

    sim.reporters.append(DCDReporter(os.path.join(OUT, f"{tag}_traj.dcd"),
                                     int(REPORT_PS * 1000 / DT)))
    per = int(REPORT_PS * 1000 / DT)
    nfr = int(ns * 1e6 / DT) // per

    rows = []
    for f in range(nfr):
        sim.step(per)
        c = np.array(sim.context.getState(getPositions=True)
                     .getPositions().value_in_unit(nanometer))
        irmsd = rmsd_after_fit(c[lig_ca], ref_lig_ca, c[iface_rec_ca], ref_if_ca) * 10
        brmsd = kabsch_rmsd(c[lig_ca], ref_lig_ca) * 10
        dm = np.linalg.norm(c[rec_heavy][:, None, :] - c[lig_heavy2][None, :, :], axis=-1)
        b = c[lig_ca] - c[lig_ca].mean(0)
        rows.append({"time_ps": (f + 1) * REPORT_PS, "iRMSD_A": irmsd,
                     "binder_RMSD_A": brmsd, "contacts": int((dm < CONTACT_CUT).sum()),
                     "min_dist_A": float(dm.min()) * 10,
                     "binder_Rg_A": float(np.sqrt((b ** 2).sum(1).mean())) * 10})
        if (f + 1) % 20 == 0 or f == nfr - 1:
            r = rows[-1]
            say(f"  {r['time_ps']/1000:.2f} ns | iRMSD {r['iRMSD_A']:.2f} | "
                f"contacts {r['contacts']} | mind {r['min_dist_A']:.2f} | "
                f"{(time.time()-t0)/60:.1f} min")

    with open(os.path.join(OUT, f"{tag}_md_metrics.json"), "w") as fh:
        json.dump({"design": tag, "ns": ns, "temp_K": TEMP, "timestep_fs": DT,
                   "solvent": "GBn2 implicit", "n_atoms": len(atoms),
                   "frames": rows}, fh, indent=2)
    with open(os.path.join(OUT, f"{tag}_final.pdb"), "w") as fh:
        PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(), fh)
    say(f"DONE {tag} in {(time.time()-t0)/60:.1f} min")
    log.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
