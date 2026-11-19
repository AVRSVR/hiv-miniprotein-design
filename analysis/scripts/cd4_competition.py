"""
Does each designed binder block CD4?

Each design was co-folded onto gp120. 1GC1 is the canonical gp120(core)-CD4(D1D2)-17b
ternary complex in the same HXB2 numbering, so superposing our gp120 onto 1GC1 gp120
places CD4 (and the CD4-induced 17b epitope) into each design's frame. Steric overlap
between the design and CD4 then reports CD4 competition directly - the neutralization
mechanism - without any docking.

Also tests pairwise co-binding: can two designs occupy gp120 simultaneously?
"""
import os, json, glob, itertools, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, Select

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SDIR = os.path.join(PROJECT, "paper_analysis", "01_structures")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")
os.makedirs(OUT, exist_ok=True)

REF = os.path.join(SDIR, "1GC1.pdb")
CLASH_HARD = 3.0     # A - severe steric overlap
CLASH_CONTACT = 4.5  # A - occupies the same surface

parser = PDBParser(QUIET=True)


def ca_map(chain):
    return {r.id[1]: r["CA"].coord for r in chain
            if r.id[0] == " " and "CA" in r}


def superpose(mob_pts, ref_pts):
    """Return (R, t) mapping mob -> ref."""
    mc, rc = mob_pts.mean(0), ref_pts.mean(0)
    m, r = mob_pts - mc, ref_pts - rc
    V, S, Wt = np.linalg.svd(m.T @ r)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, rc - (mc @ R)


def heavy_coords(entity):
    return np.array([a.coord for a in entity.get_atoms()
                     if a.element != "H"])


ref_st = parser.get_structure("ref", REF)
ref_model = next(iter(ref_st))
ref_gp120 = ref_model["G"]
ref_cd4 = ref_model["C"]
ref_gp120_ca = ca_map(ref_gp120)
cd4_xyz_ref = heavy_coords(ref_cd4)
# 17b Fab (CD4-induced epitope antibody)
b17_xyz_ref = np.vstack([heavy_coords(ref_model[c]) for c in ("L", "H")
                         if c in ref_model])

print(f"1GC1 reference: gp120 {len(ref_gp120_ca)} CA, CD4 {len(cd4_xyz_ref)} heavy atoms")

# gp120 residues forming the CD4 epitope in 1GC1 (any heavy atom within 5 A of CD4)
cd4_epitope_res = set()
for r in ref_gp120:
    if r.id[0] != " ":
        continue
    ra = np.array([a.coord for a in r if a.element != "H"])
    if ra.size and np.linalg.norm(ra[:, None, :] - cd4_xyz_ref[None, :, :],
                                  axis=-1).min() < 5.0:
        cd4_epitope_res.add(r.id[1])
print(f"CD4 epitope on gp120: {len(cd4_epitope_res)} residues "
      f"({sorted(cd4_epitope_res)[:12]}...)")

rows = []
placed_binders = {}   # design -> binder heavy coords in the 1GC1 frame

for f in sorted(glob.glob(os.path.join(SDIR, "*_split.pdb"))):
    tag = os.path.basename(f).replace("_split.pdb", "")
    st = parser.get_structure(tag, f)
    m = next(iter(st))
    if "A" not in m or "B" not in m:
        continue
    gp = ca_map(m["A"])

    common = sorted(set(gp) & set(ref_gp120_ca))
    if len(common) < 30:
        print(f"[skip] {tag}: only {len(common)} shared gp120 residues")
        continue

    # Fit locally on the CD4 epitope rather than all of gp120: our complexes are
    # full-length gp120 while 1GC1 is the core, so a global fit carries loop
    # differences into the region the competition question actually depends on.
    local = [i for i in common if i in cd4_epitope_res]
    use = local if len(local) >= 15 else common
    mob = np.array([gp[i] for i in use])
    ref = np.array([ref_gp120_ca[i] for i in use])
    R, t = superpose(mob, ref)
    fit_rmsd = float(np.sqrt((((mob @ R) + t - ref) ** 2).sum(1).mean()))

    glob_mob = np.array([gp[i] for i in common])
    glob_ref = np.array([ref_gp120_ca[i] for i in common])
    global_rmsd = float(np.sqrt((((glob_mob @ R) + t - glob_ref) ** 2).sum(1).mean()))

    binder = heavy_coords(m["B"])
    binder_ref_frame = (binder @ R) + t
    placed_binders[tag] = binder_ref_frame

    d_cd4 = np.linalg.norm(binder_ref_frame[:, None, :] - cd4_xyz_ref[None, :, :], axis=-1)
    hard = int((d_cd4 < CLASH_HARD).sum())
    cont = int((d_cd4 < CLASH_CONTACT).sum())
    mind = float(d_cd4.min())

    d_17b = np.linalg.norm(binder_ref_frame[:, None, :] - b17_xyz_ref[None, :, :], axis=-1)
    hard17 = int((d_17b < CLASH_HARD).sum())
    mind17 = float(d_17b.min())

    verdict = ("CD4-blocking" if hard >= 20 else
               "partial overlap" if cont >= 20 else
               "non-competitive")

    rows.append({
        "design": tag,
        "gp120_superposition_residues": len(common),
        "cd4_epitope_fit_RMSD_A": round(fit_rmsd, 2),
        "global_gp120_RMSD_A": round(global_rmsd, 2),
        "binder_CD4_min_dist_A": round(mind, 2),
        "binder_CD4_clash_atoms_3A": hard,
        "binder_CD4_contact_atoms_4.5A": cont,
        "binder_17b_min_dist_A": round(mind17, 2),
        "binder_17b_clash_atoms_3A": hard17,
        "CD4_competition": verdict,
    })
    print(f"[ok] {tag:<16} localfit {fit_rmsd:5.2f} (glob {global_rmsd:5.2f}) | "
          f"CD4 min {mind:6.2f} A | clash {hard:5d} | contact {cont:5d} -> {verdict}")

# pairwise co-binding: two designs on gp120 at once
pairs = []
for a, b in itertools.combinations(sorted(placed_binders), 2):
    d = np.linalg.norm(placed_binders[a][:, None, :] - placed_binders[b][None, :, :], axis=-1)
    hard = int((d < CLASH_HARD).sum())
    pairs.append({"design_1": a, "design_2": b,
                  "min_dist_A": round(float(d.min()), 2),
                  "clash_atoms_3A": hard,
                  "co_binding": "compatible" if hard == 0 else "mutually exclusive"})

with open(os.path.join(OUT, "cd4_competition.json"), "w") as fh:
    json.dump({"reference": "1GC1 (gp120 core + CD4 D1D2 + 17b Fab)",
               "clash_cutoffs_A": {"hard": CLASH_HARD, "contact": CLASH_CONTACT},
               "per_design": rows, "pairwise_cobinding": pairs}, fh, indent=2)

print("\n=== pairwise co-binding on gp120 ===")
for p in sorted(pairs, key=lambda x: x["clash_atoms_3A"]):
    print(f"{p['design_1']:<16} + {p['design_2']:<16} min {p['min_dist_A']:6.2f} A  "
          f"clash {p['clash_atoms_3A']:5d}  {p['co_binding']}")
