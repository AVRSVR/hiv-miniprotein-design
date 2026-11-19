"""
Does each binder's footprint land on surface that an N-glycan would cover?

Every structure used here has its glycans stripped: heteroatoms were removed
during preparation, so all interfaces in this work are computed on a naked gp120
surface. gp120 is in reality one of the most heavily glycosylated proteins known,
and a designed binder that lands on glycan-shielded surface would be unable to
reach its epitope no matter how good the protein-protein energetics look.

That risk is bounded rather than assumed here. N-linked sequons (N-X-S/T, X != P)
are found from the model's own sequence and numbering, requiring consecutive
residue numbers so detection never jumps a chain break. For each footprint we
report the nearest sequon asparagine and how many sequons fall within 10 and 15 A
- roughly the reach of the core pentasaccharide and of a processed complex glycan.

The controls are the point: CD4 and the 17b Fab bind fully glycosylated gp120 in
the virion, so their footprints measure how much glycan proximity a working
partner actually tolerates. A design inside that range is not excused of the
limitation, but it is no worse placed than the natural ligands.

Writes glycan_occlusion.json and table14_glycan_occlusion.csv.
"""
import os, csv, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCENE = os.path.join(PA, "04_figures", "scene")
OUT = os.path.join(PA, "05_tables")
CUT = 5.0            # A, footprint definition, as in footprint_overlap.py
NEAR = (10.0, 15.0)  # A, core-glycan and complex-glycan reach

NAME = {"S310": "D1 (CD4-mimetic control)", "Chai1_Best": "D2 (de novo bispecific)",
        "Cyclic_Ring": "D3 (CD4-mimetic cyclic)", "model_4_chai": "D2' (alternative pose)",
        "Latest_Best": "D5 (de novo, coreceptor site)"}

p = PDBParser(QUIET=True)


def residues(model, chains=None):
    out = []
    for ch in model:
        if chains and ch.id not in chains:
            continue
        for r in ch:
            if r.id[0] == " ":
                out.append(r)
    return out


def coords(res_list):
    return [(f"{r.get_parent().id}{r.id[1]}", a.coord)
            for r in res_list for a in r]


def footprint(target_res, probe_res):
    P = np.array([c for _, c in coords(probe_res)])
    keep = []
    for r in target_res:
        X = np.array([a.coord for a in r])
        if np.min(np.linalg.norm(X[:, None, :] - P[None, :, :], axis=-1)) < CUT:
            keep.append(r)
    return keep


gp = residues(p.get_structure("g", os.path.join(SCENE, "gp120_ref.pdb"))[0])
by_num = {r.id[1]: r for r in gp}

# ---- N-linked sequons from the model's own sequence, never crossing a gap
THREE1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}
sequons = []
for n in sorted(by_num):
    r0, r1, r2 = by_num.get(n), by_num.get(n + 1), by_num.get(n + 2)
    if not (r0 and r1 and r2):
        continue
    a0, a1, a2 = (THREE1.get(x.get_resname(), "X") for x in (r0, r1, r2))
    if a0 == "N" and a1 != "P" and a2 in ("S", "T"):
        sequons.append(n)

print("gp120 model: %d residues, numbering %d-%d" % (len(gp), min(by_num), max(by_num)))
print("N-linked sequons found (%d): %s" % (len(sequons), ", ".join("N%d" % n for n in sequons)))

seq_atoms = {n: np.array([a.coord for a in by_num[n]]) for n in sequons}


def profile(label, fp_res):
    """Nearest sequon and sequon counts within each reach, for one footprint."""
    if not fp_res:
        return None
    F = np.array([a.coord for r in fp_res for a in r])
    ids = {r.id[1] for r in fp_res}
    d = {n: float(np.min(np.linalg.norm(F[:, None, :] - X[None, :, :], axis=-1)))
         for n, X in seq_atoms.items()}
    on = sorted(ids & set(sequons))
    nearest = min(d, key=d.get) if d else None
    row = {
        "label": label,
        "footprint_res": len(fp_res),
        "sequon_asn_in_footprint": ["N%d" % n for n in on],
        "nearest_sequon": "N%d" % nearest if nearest else "",
        "nearest_sequon_dist_A": round(d[nearest], 2) if nearest else None,
    }
    for cut in NEAR:
        within = sorted(n for n, v in d.items() if v < cut)
        row["sequons_within_%gA" % cut] = len(within)
        row["sequons_within_%gA_list" % cut] = ["N%d" % n for n in within]
    return row


ref = p.get_structure("r", os.path.join(PA, "01_structures", "1GC1.pdb"))[0]
rows = []
rows.append(profile("CD4 (experimental, binds glycosylated gp120)",
                    footprint(gp, residues(ref, {"C"}))))
rows.append(profile("17b Fab (experimental, binds glycosylated gp120)",
                    footprint(gp, residues(ref, {"L", "H"}))))
for tag, nm in NAME.items():
    f = os.path.join(SCENE, "binder_%s.pdb" % tag)
    if os.path.exists(f):
        rows.append(profile(nm, footprint(gp, residues(p.get_structure(tag, f)[0]))))

rows = [r for r in rows if r]
print()
hdr = "%-46s %4s %6s %8s %7s %7s" % ("footprint", "res", "on-seq", "nearest", "<10A", "<15A")
print(hdr)
print("-" * len(hdr))
for r in rows:
    print("%-46s %4d %6d %8s %7d %7d" % (
        r["label"][:46], r["footprint_res"], len(r["sequon_asn_in_footprint"]),
        "%.1f" % r["nearest_sequon_dist_A"], r["sequons_within_10A"],
        r["sequons_within_15A"]))

with open(os.path.join(OUT, "glycan_occlusion.json"), "w") as fh:
    json.dump({"footprint_cutoff_A": CUT, "reach_A": list(NEAR),
               "sequons": ["N%d" % n for n in sequons], "per_footprint": rows},
              fh, indent=2)

with open(os.path.join(OUT, "table14_glycan_occlusion.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Footprint", "Footprint res", "Sequon Asn in footprint",
                "Nearest sequon", "Distance (A)", "Sequons <10 A", "Sequons <15 A",
                "Sequons <15 A (list)"])
    for r in rows:
        w.writerow([r["label"], r["footprint_res"],
                    ", ".join(r["sequon_asn_in_footprint"]) or "none",
                    r["nearest_sequon"], r["nearest_sequon_dist_A"],
                    r["sequons_within_10A"], r["sequons_within_15A"],
                    ", ".join(r["sequons_within_15A_list"])])
print("\nwrote glycan_occlusion.json and table14_glycan_occlusion.csv")
