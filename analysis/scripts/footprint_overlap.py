"""
How much of each receptor's gp120 footprint does a binder actually cover?

Steric clash counts against CD4 or the 17b Fab answer "would the receptor be
blocked?", but they do not say whether the binder sits *on* the receptor's
epitope or merely collides with the bulk of the receptor from an adjacent
surface. Those are different mechanistic claims and only the second is supported
when footprint overlap is low, so the overlap is measured here explicitly.

A footprint is the set of gp120 residues within 5 A of the partner in the
superposed 1GC1 frame.
"""
import os, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

PA = r"C:\Users\anant\Desktop\hiv_miniprotein_design\paper_analysis"
SCENE = os.path.join(PA, "04_figures", "scene")
OUT = os.path.join(PA, "05_tables")
CUT = 5.0
NAME = {"S310": "D1 (CD4-mimetic control)", "Chai1_Best": "D2 (de novo bispecific)",
        "Cyclic_Ring": "D3 (CD4-mimetic cyclic)", "model_4_chai": "D2' (alternative pose)",
        "Latest_Best": "D5 (de novo, coreceptor site)"}

p = PDBParser(QUIET=True)


def atoms(model, chains=None):
    out = []
    for ch in model:
        if chains and ch.id not in chains:
            continue
        for r in ch:
            if r.id[0] != " ":
                continue
            for a in r:
                out.append((f"{ch.id}{r.id[1]}", a.coord))
    return out


def footprint(target, probe):
    P = np.array([c for _, c in probe])
    return {rid for rid, c in target
            if np.min(np.linalg.norm(P - c, axis=1)) < CUT}


gp = atoms(p.get_structure("g", os.path.join(SCENE, "gp120_ref.pdb"))[0])
ref = p.get_structure("r", os.path.join(PA, "01_structures", "1GC1.pdb"))[0]
cd4_fp = footprint(gp, atoms(ref, {"C"}))
ab_fp = footprint(gp, atoms(ref, {"L", "H"}))

rows = []
for tag in NAME:
    f = os.path.join(SCENE, f"binder_{tag}.pdb")
    if not os.path.exists(f):
        continue
    b = footprint(gp, atoms(p.get_structure(tag, f)[0]))
    rows.append({
        "design": tag, "label": NAME[tag],
        "binder_footprint_res": len(b),
        "cd4_footprint_res": len(cd4_fp), "17b_footprint_res": len(ab_fp),
        "shared_with_CD4": len(b & cd4_fp),
        "shared_with_17b": len(b & ab_fp),
        "pct_of_CD4_footprint": round(100 * len(b & cd4_fp) / len(cd4_fp), 1),
        "pct_of_17b_footprint": round(100 * len(b & ab_fp) / len(ab_fp), 1),
    })
    r = rows[-1]
    print(f"{NAME[tag]:<34} own {r['binder_footprint_res']:>3} res | "
          f"CD4 {r['pct_of_CD4_footprint']:>5.1f}% | 17b {r['pct_of_17b_footprint']:>5.1f}%")

with open(os.path.join(OUT, "footprint_overlap.json"), "w") as fh:
    json.dump({"cutoff_A": CUT, "per_design": rows}, fh, indent=2)
print("\nwrote footprint_overlap.json")
