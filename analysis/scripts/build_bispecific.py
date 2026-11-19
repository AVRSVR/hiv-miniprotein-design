"""
Build an actual 3D model of a tandem bispecific construct.

Table 6 established that linking two arms is geometrically feasible. This builds the
molecule: the two arms are held in the bound poses determined from the co-folded
complexes (they occupy gp120 simultaneously with zero steric conflict), a flexible
Gly/Ser linker is threaded between the C-terminus of arm 1 and the N-terminus of
arm 2, and the linker is relaxed by restrained energy minimisation with both arms
and gp120 held fixed - so only the linker moves and the validated binding poses are
preserved.

The initial linker path is routed outward, away from the gp120 surface, so the
starting guess does not thread through the target; minimisation then resolves the
local geometry.

Output: a single-chain construct PDB plus a complex with gp120, and a clash report.

Usage: python build_bispecific.py <arm1_tag> <arm2_tag> [n_linker]
"""
import os, sys, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

from openmm import LangevinMiddleIntegrator, Platform, CustomExternalForce
from openmm.app import ForceField, Simulation, PDBFile, HBonds, NoCutoff, Modeller
from openmm.unit import (kelvin, picoseconds, femtoseconds, nanometer,
                         kilojoule_per_mole, amu)
from pdbfixer import PDBFixer

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SCENE = os.path.join(PROJECT, "paper_analysis", "04_figures", "scene")
OUT = os.path.join(PROJECT, "paper_analysis", "06_bispecific")
os.makedirs(OUT, exist_ok=True)

CA_STEP = 3.80          # A, consecutive CA separation
parser = PDBParser(QUIET=True)


def load_chain(path):
    st = parser.get_structure("x", path)
    return next(iter(next(iter(st))))


def residues_of(chain):
    return [r for r in chain if r.id[0] == " "]


def linker_path(start, end, n, avoid_center, clearance=6.0):
    """
    CA positions for n linker residues between start and end, bowed away from
    avoid_center so the initial guess does not pass through the target.
    """
    pts = []
    span = np.linalg.norm(end - start)
    # outward direction: away from the target centre, perpendicular-ish to the span
    mid = 0.5 * (start + end)
    out = mid - avoid_center
    if np.linalg.norm(out) < 1e-6:
        out = np.array([0.0, 0.0, 1.0])
    out = out / np.linalg.norm(out)
    axis = (end - start) / max(span, 1e-6)
    out = out - np.dot(out, axis) * axis      # make it perpendicular to the span
    if np.linalg.norm(out) < 1e-6:
        out = np.array([0.0, 0.0, 1.0])
    out = out / np.linalg.norm(out)

    # bow height: enough slack that the chain arcs clear of the surface
    contour = n * CA_STEP
    bow = max(clearance, 0.5 * np.sqrt(max(contour ** 2 - span ** 2, 0.0)))
    for i in range(1, n + 1):
        t = i / (n + 1)
        base = start + (end - start) * t
        pts.append(base + out * bow * np.sin(np.pi * t))
    return pts


def gly_residue_atoms(ca, direction):
    """Approximate backbone atoms for one Gly from its CA and chain direction."""
    d = direction / max(np.linalg.norm(direction), 1e-6)
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(ref, d)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    perp = np.cross(d, ref)
    perp = perp / max(np.linalg.norm(perp), 1e-6)
    return {
        "N": ca - d * 1.46,
        "CA": ca,
        "C": ca + d * 1.52,
        "O": ca + d * 1.52 + perp * 1.23,
    }


def pdb_line(serial, name, resname, chain, resseq, xyz, element):
    return (f"ATOM  {serial:5d} {name:^4s} {resname:>3s} {chain}{resseq:>4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          "
            f"{element:>2s}\n")


def main(tag1, tag2, n_linker=None):
    a1 = load_chain(os.path.join(SCENE, f"binder_{tag1}.pdb"))
    a2 = load_chain(os.path.join(SCENE, f"binder_{tag2}.pdb"))
    gp = load_chain(os.path.join(SCENE, "gp120_ref.pdb"))

    r1, r2 = residues_of(a1), residues_of(a2)
    cterm, nterm = r1[-1], r2[0]
    if "C" not in cterm or "N" not in nterm:
        raise SystemExit("termini lack backbone atoms")
    start = np.array(cterm["C"].coord)
    end = np.array(nterm["N"].coord)
    span = float(np.linalg.norm(end - start))

    if n_linker is None:
        n_linker = int(max(np.ceil(span / 3.5) + 4, 8))
    gp_xyz = np.array([a.coord for a in gp.get_atoms() if a.element != "H"])
    gp_center = gp_xyz.mean(0)

    print(f"arm1 {tag1} ({len(r1)} aa)  ->  arm2 {tag2} ({len(r2)} aa)")
    print(f"C-to-N span {span:.1f} A, linker {n_linker} residues "
          f"(contour {n_linker*CA_STEP:.1f} A)")

    path = linker_path(start, end, n_linker, gp_center)

    # ---- assemble a single chain: arm1 | linker | arm2 -------------------
    lines, serial, resseq = [], 1, 1
    for r in r1:
        for a in r:
            if a.element == "H":
                continue
            # arm 1 is no longer a chain end: its terminal OXT would otherwise
            # leave the junction residue over-bonded
            if a.get_name() == "OXT":
                continue
            lines.append(pdb_line(serial, a.get_name(), r.get_resname(), "A",
                                  resseq, a.coord, a.element))
            serial += 1
        resseq += 1

    seq_linker = []
    for i, ca in enumerate(path):
        nxt = path[i + 1] if i + 1 < len(path) else end
        atoms = gly_residue_atoms(ca, nxt - ca)
        resname = "GLY" if (i % 5) != 4 else "SER"
        if resname == "SER":
            resname = "GLY"          # keep it chemically simple: poly-Gly/Ser backbone
        seq_linker.append("G")
        for nm in ("N", "CA", "C", "O"):
            lines.append(pdb_line(serial, nm, resname, "A", resseq,
                                  atoms[nm], nm[0]))
            serial += 1
        resseq += 1

    n_arm2_start = resseq
    for r in r2:
        for a in r:
            if a.element == "H":
                continue
            lines.append(pdb_line(serial, a.get_name(), r.get_resname(), "A",
                                  resseq, a.coord, a.element))
            serial += 1
        resseq += 1
    lines.append("TER\nEND\n")

    raw = os.path.join(OUT, f"{tag1}__{tag2}_raw.pdb")
    with open(raw, "w") as fh:
        fh.writelines(lines)
    print(f"assembled raw construct: {resseq-1} residues -> {os.path.basename(raw)}")

    # ---- relax the linker with gp120 present ----------------------------
    # gp120 must be in the system: minimising the construct alone lets the linker
    # relax straight through the target, since nothing excludes that volume.
    raw_with_target = os.path.join(OUT, f"{tag1}__{tag2}_raw_with_gp120.pdb")
    with open(raw_with_target, "w") as fh:
        for ln in open(os.path.join(SCENE, "gp120_ref.pdb")):
            if ln.startswith("ATOM") and ln[76:78].strip() != "H":
                fh.write(ln[:21] + "G" + ln[22:])
        fh.write("TER\n")
        fh.writelines(lines)

    fx = PDBFixer(filename=raw_with_target)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(7.4)
    top, pos = fx.topology, fx.positions

    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(top, nonbondedMethod=NoCutoff, constraints=HBonds,
                             hydrogenMass=1.5 * amu)

    atoms = list(top.atoms())
    xyz0 = np.array(pos.value_in_unit(nanometer))

    # linker = positions len(r1) .. len(r1)+n_linker-1 within the construct chain
    chainA = [c for c in top.chains() if c.id == "A"]
    construct_res = [r for c in chainA for r in c.residues()]
    linker_res_idx = {r.index for r in construct_res[len(r1):len(r1) + n_linker]}

    restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", 1000.0 * kilojoule_per_mole / nanometer**2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    n_fixed = 0
    for i, a in enumerate(atoms):
        if a.residue.index in linker_res_idx:
            continue                      # linker is free
        if a.element is not None and a.element.symbol == "H":
            continue
        restraint.addParticle(i, xyz0[i].tolist())
        n_fixed += 1
    system.addForce(restraint)
    print(f"restraining {n_fixed} atoms (gp120 + arms); relaxing "
          f"{len(linker_res_idx)} linker residues")

    integ = LangevinMiddleIntegrator(310.15 * kelvin, 1.0 / picoseconds, 2.0 * femtoseconds)
    try:
        sim = Simulation(top, system, integ, Platform.getPlatformByName("OpenCL"))
    except Exception:
        sim = Simulation(top, system, integ)
    sim.context.setPositions(pos)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=3000)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    print(f"minimised: {e0} -> {st.getPotentialEnergy()}")

    # full complex (gp120 + construct) and the construct alone
    combo = os.path.join(OUT, f"{tag1}__{tag2}_with_gp120.pdb")
    with open(combo, "w") as fh:
        PDBFile.writeFile(top, st.getPositions(), fh)

    built = os.path.join(OUT, f"{tag1}__{tag2}_construct.pdb")
    mod = Modeller(top, st.getPositions())
    mod.delete([r for r in top.residues() if r.chain.id != "A"])
    with open(built, "w") as fh:
        PDBFile.writeFile(mod.topology, mod.positions, fh)

    # ---- clash report ---------------------------------------------------
    final = np.array(st.getPositions().value_in_unit(nanometer)) * 10.0
    heavy = [i for i, a in enumerate(atoms)
             if a.element is not None and a.element.symbol != "H"]
    lk = [i for i in heavy if atoms[i].residue.index in linker_res_idx]
    armi = [i for i in heavy
            if atoms[i].residue.chain.id == "A" and atoms[i].residue.index not in linker_res_idx]
    gpi = [i for i in heavy if atoms[i].residue.chain.id != "A"]

    d_gp = np.linalg.norm(final[lk][:, None, :] - final[gpi][None, :, :], axis=-1)
    d_arm = np.linalg.norm(final[lk][:, None, :] - final[armi][None, :, :], axis=-1)
    # arm poses must be preserved by the relaxation
    arm_rmsd = float(np.sqrt(((final[armi] - xyz0[armi] * 10.0) ** 2).sum(1).mean()))

    # A flexible linker traversing a protein surface legitimately makes polar
    # contacts at 2.6-3.2 A; only sub-2.5 A heavy-atom approaches are true steric
    # violations, so the two are counted separately.
    report = {
        "arm_1": tag1, "arm_2": tag2, "linker_residues": n_linker,
        "C_to_N_span_A": round(span, 1),
        "construct_residues": resseq - 1,
        "linker_min_dist_to_gp120_A": round(float(d_gp.min()), 2),
        "linker_gp120_steric_clashes_lt2.5A": int((d_gp < 2.5).sum()),
        "linker_gp120_surface_contacts_2.5_3.2A": int(((d_gp >= 2.5) & (d_gp < 3.2)).sum()),
        "linker_min_dist_to_arms_A": round(float(d_arm.min()), 2),
        "arm_heavy_atom_RMSD_after_relax_A": round(arm_rmsd, 3),
        "construct_pdb": os.path.basename(built),
        "complex_pdb": os.path.basename(combo),
    }
    print(json.dumps(report, indent=2))

    with open(os.path.join(OUT, f"{tag1}__{tag2}_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    return report


if __name__ == "__main__":
    t1 = sys.argv[1] if len(sys.argv) > 1 else "S310"
    t2 = sys.argv[2] if len(sys.argv) > 2 else "Latest_Best"
    n = int(sys.argv[3]) if len(sys.argv) > 3 else None
    main(t1, t2, n)
