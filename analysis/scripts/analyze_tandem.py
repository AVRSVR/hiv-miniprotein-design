"""
Per-arm analysis of the tandem construct trajectory.

run_md2.py reports metrics for the binder chain as a whole, which for a two-armed
construct cannot answer the question the run was launched to settle: does each arm
stay on its own epitope once the linker is free and nothing is restrained? A whole
construct iRMSD can look excellent while one arm has let go and the other has
rotated to compensate.

The construct is chain B, 185 residues in N-to-C order: arm 1 is D1 (97 residues),
then a 22-residue poly-glycine linker, then arm 2 is D5 (66 residues). Each arm is
scored separately against the receptor, each superposed on the receptor interface
residues nearest to it so that one arm's drift cannot contaminate the other's iRMSD.

Also checked: whether the linker stays clear of gp120, and whether the compressed
CA-CA distances the builder left in it (four near 3.1 A against a trans value of
3.8 A) survive dynamics - which says whether the published construct model should be
replaced by the relaxed one.

Run: python analyze_tandem.py
"""
import io
import json
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
import mdtraj as md
from pdbfixer import PDBFixer
from openmm.app import PDBFile

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(PA, "03_md")
OUT = os.path.join(PA, "05_tables")
TAG = "Tandem_D1_D5"

ARM1, LINK, ARM2 = 97, 22, 66          # residues, N to C
CONTACT = 0.45                          # nm
NEAR = 0.8                              # nm, receptor residues used for superposition
TOPO = os.path.join(MD, f"{TAG}_topology.pdb")


def topology():
    """Rebuild the simulated topology deterministically from the truncated PDB.

    The DCD was written for the PDBFixer-completed topology, not for the truncated
    file on disk, so the same fixer calls are repeated here. The truncated file was
    already written after short fragments were dropped, so this reproduces the
    simulated atom set exactly.
    """
    if not os.path.exists(TOPO):
        fx = PDBFixer(filename=os.path.join(MD, f"{TAG}_trunc.pdb"))
        fx.findMissingResidues()
        fx.missingResidues = {}
        fx.findMissingAtoms()
        fx.addMissingAtoms()
        fx.addMissingHydrogens(7.4)
        with io.open(TOPO, "w") as fh:
            PDBFile.writeFile(fx.topology, fx.positions, fh)
    return md.load(TOPO)


ref = topology()
t = md.load(os.path.join(MD, f"{TAG}_traj.dcd"), top=ref.topology)
print(f"{t.n_frames} frames, {t.n_atoms} atoms")
assert t.n_atoms == ref.n_atoms, (t.n_atoms, ref.n_atoms)

top = t.topology
binder_res = [r for r in top.residues if r.chain.index == next(
    c.index for c in top.chains
    if any(a.residue.name != "HOH" for a in c.atoms)
    and sum(1 for _ in c.residues) == ARM1 + LINK + ARM2)]
assert len(binder_res) == ARM1 + LINK + ARM2, len(binder_res)

groups = {
    "arm1_D1": binder_res[:ARM1],
    "linker": binder_res[ARM1:ARM1 + LINK],
    "arm2_D5": binder_res[ARM1 + LINK:],
}
rec_res = [r for r in top.residues if r not in set(binder_res)]


def heavy(res_list):
    return [a.index for r in res_list for a in r.atoms if a.element.symbol != "H"]


def ca(res_list):
    return [a.index for r in res_list for a in r.atoms if a.name == "CA"]


rec_heavy = heavy(rec_res)
X0 = t.xyz[0]

rows = {}
for name, res in groups.items():
    g_heavy = heavy(res)
    # receptor residues starting near this group, used as the superposition frame
    d0 = np.linalg.norm(X0[rec_heavy][:, None, :] - X0[g_heavy][None, :, :], axis=-1)
    near_idx = {top.atom(rec_heavy[i]).residue.index
                for i in np.where(d0.min(axis=1) < NEAR)[0]}
    near_ca = [a.index for r in rec_res if r.index in near_idx
               for a in r.atoms if a.name == "CA"]

    contacts, mind = [], []
    for f in range(t.n_frames):
        D = np.linalg.norm(t.xyz[f][rec_heavy][:, None, :]
                           - t.xyz[f][g_heavy][None, :, :], axis=-1)
        contacts.append(int((D < CONTACT).sum()))
        mind.append(float(D.min()))

    row = {"n_res": len(res), "contacts_start": contacts[0],
           "contacts_mean": float(np.mean(contacts)),
           "contacts_final_half": float(np.mean(contacts[len(contacts) // 2:])),
           "contacts_retained_pct": round(100.0 * np.mean(contacts[len(contacts) // 2:])
                                          / max(contacts[0], 1), 1),
           "min_dist_mean_A": round(10 * float(np.mean(mind)), 2),
           "min_dist_final_A": round(10 * mind[-1], 2),
           "bound_at_end": bool(mind[-1] < 0.5),
           "n_near_receptor_ca": len(near_ca)}

    if len(near_ca) >= 4 and name != "linker":
        sel = np.array(near_ca + ca(res))
        sub = t.atom_slice(sel)
        sub.superpose(sub, 0, atom_indices=np.arange(len(near_ca)))
        g = np.arange(len(near_ca), sub.n_atoms)
        dev = np.sqrt(((sub.xyz[:, g, :] - sub.xyz[0, g, :]) ** 2).sum(-1).mean(-1))
        row["iRMSD_mean_A"] = round(10 * float(dev.mean()), 2)
        row["iRMSD_final_half_A"] = round(10 * float(dev[len(dev) // 2:].mean()), 2)
        row["iRMSD_max_A"] = round(10 * float(dev.max()), 2)
        fold = md.rmsd(t.atom_slice(np.array(ca(res))),
                       t.atom_slice(np.array(ca(res))), 0)
        row["fold_RMSD_final_half_A"] = round(10 * float(np.mean(fold[len(fold) // 2:])), 2)
    rows[name] = row

# did the builder's compressed peptide geometry relax?
link_ca = ca(groups["linker"])
def ca_ca(frame):
    P = frame[link_ca]
    return np.linalg.norm(P[1:] - P[:-1], axis=1) * 10.0
first, last = ca_ca(t.xyz[0]), ca_ca(t.xyz[-1])
geom = {"linker_CA_CA_start_min_A": round(float(first.min()), 2),
        "linker_CA_CA_start_below_3.4A": int((first < 3.4).sum()),
        "linker_CA_CA_end_min_A": round(float(last.min()), 2),
        "linker_CA_CA_end_below_3.4A": int((last < 3.4).sum()),
        "linker_CA_CA_end_mean_A": round(float(last.mean()), 2)}

print()
hdr = "%-10s %4s %9s %9s %8s %8s %8s %7s" % (
    "group", "res", "contacts0", "contactsF", "retain%", "mind(A)", "iRMSD_F", "bound")
print(hdr)
print("-" * len(hdr))
for name, r in rows.items():
    print("%-10s %4d %9d %9.0f %8.1f %8.2f %8s %7s" % (
        name, r["n_res"], r["contacts_start"], r["contacts_final_half"],
        r["contacts_retained_pct"], r["min_dist_final_A"],
        r.get("iRMSD_final_half_A", "-"), r["bound_at_end"]))
print("\nlinker geometry:", json.dumps(geom))

with io.open(os.path.join(OUT, "tandem_md.json"), "w", encoding="utf-8") as fh:
    json.dump({"tag": TAG, "n_frames": t.n_frames,
               "protocol": "2 fs timestep, no hydrogen mass repartitioning, "
                           "20000-step minimisation, 4x20 ps heating ramp",
               "per_group": rows, "linker_geometry": geom}, fh, indent=2)
print("wrote tandem_md.json")
