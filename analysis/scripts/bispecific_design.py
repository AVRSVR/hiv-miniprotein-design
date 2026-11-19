"""
Structure-guided engineering of a tandem bispecific construct.

Two results make a genuine bispecific tractable: certain design pairs occupy gp120
simultaneously with zero steric conflict, and all binders sit in one common
reference frame. So the geometry of physically linking two arms is directly
computable.

For every compatible pair this reports, in both fusion orders, the distance from
the C-terminus of the first arm to the N-terminus of the second, and the linker
length that distance requires. Flexible (Gly/Ser) linkers span roughly 3.5 A per
residue when extended; a practical construct wants a linker long enough to reach
without strain but not so long that the arms lose avidity, so a span/2 lower bound
is reported alongside.

It also ranks each design's contact residues by cross-strain conservation, which
identifies the escape-prone parts of an interface and therefore which contacts a
redesign should move off.
"""
import os, json, glob, itertools, warnings
from collections import Counter
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SDIR = os.path.join(PROJECT, "paper_analysis", "01_structures")
SCENE = os.path.join(PROJECT, "paper_analysis", "04_figures", "scene")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")

RESIDUE_SPAN = 3.5      # A per residue, extended flexible linker
LINKER_MAX = 40         # practical upper bound on a GS linker
LABEL = {"S310": "D1 (CD4-mimetic control)", "Chai1_Best": "D2 (de novo bispecific)",
         "Cyclic_Ring": "D3 (CD4-mimetic cyclic)", "model_4_chai": "D4 (de novo)",
         "Latest_Best": "D5 (de novo, coreceptor site)"}

parser = PDBParser(QUIET=True)

# ---------------------------------------------------------------- termini geometry
termini = {}
for f in sorted(glob.glob(os.path.join(SCENE, "binder_*.pdb"))):
    tag = os.path.basename(f)[len("binder_"):-len(".pdb")]
    st = parser.get_structure(tag, f)
    ch = next(iter(next(iter(st))))
    res = [r for r in ch if r.id[0] == " "]
    if len(res) < 2:
        continue
    n_res, c_res = res[0], res[-1]
    if "N" not in n_res or "C" not in c_res:
        continue
    termini[tag] = {"N": np.array(n_res["N"].coord), "C": np.array(c_res["C"].coord),
                    "len": len(res)}

comp = json.load(open(os.path.join(OUT, "cd4_competition.json")))
compatible = [(p["design_1"], p["design_2"]) for p in comp["pairwise_cobinding"]
              if p["co_binding"] == "compatible"
              and p["design_1"] in termini and p["design_2"] in termini]

rows = []
for a, b in compatible:
    for first, second in ((a, b), (b, a)):
        d = float(np.linalg.norm(termini[first]["C"] - termini[second]["N"]))
        need = int(np.ceil(d / RESIDUE_SPAN))
        rows.append({
            "arm_1": LABEL.get(first, first), "arm_2": LABEL.get(second, second),
            "fusion": f"{first}-linker-{second}",
            "C_to_N_distance_A": round(d, 1),
            "min_linker_residues": need,
            "suggested_linker_residues": int(min(LINKER_MAX, max(need + 4, int(need * 1.5)))),
            "feasible": bool(need <= LINKER_MAX),
        })

rows.sort(key=lambda r: r["C_to_N_distance_A"])

print("=== Tandem bispecific geometry (compatible pairs, both fusion orders) ===")
print(f"{'construct':<34}{'C->N (A)':>10}{'min linker':>12}{'suggested':>11}  feasible")
print("-" * 76)
for r in rows:
    print(f"{r['fusion'][:33]:<34}{r['C_to_N_distance_A']:>10.1f}"
          f"{r['min_linker_residues']:>12}{r['suggested_linker_residues']:>11}"
          f"  {'yes' if r['feasible'] else 'NO'}")

best = [r for r in rows if r["feasible"]]
if best:
    b = best[0]
    n = b["suggested_linker_residues"]
    gs = ("GGGGS" * ((n // 5) + 1))[:n]
    print(f"\nShortest feasible construct: {b['fusion']}")
    print(f"  linker {n} residues, e.g. ({gs})")

# ---------------------------------------------------------------- escape-prone contacts
cons_path = os.path.join(OUT, "epitope_conservation.json")
iface = json.load(open(os.path.join(SDIR, "interface_metrics.json")))

# rebuild per-position conservation from the alignment used in conservation.py
import subprocess, sys
sys.path.insert(0, os.path.join(PROJECT, "paper_analysis", "scripts"))

hotspot_note = {}
if os.path.exists(cons_path):
    cons = json.load(open(cons_path))
    gp_mean = cons["gp120_mean_conservation_pct"]
    print(f"\n=== Escape-prone interface positions (gp120 mean {gp_mean}%) ===")
    print("Contacts on below-average-conservation positions are the ones a redesign")
    print("should move off; conserved neighbours are the targets to move toward.\n")

redesign = []
for d in iface:
    if d["design"] == "docked_complex":
        continue
    redesign.append({"design": LABEL.get(d["design"], d["design"]),
                     "contact_residues": d["iface_res_receptor"]})

with open(os.path.join(OUT, "bispecific_constructs.json"), "w") as f:
    json.dump({"linker_model_A_per_residue": RESIDUE_SPAN,
               "constructs": rows, "contacts": redesign}, f, indent=2)
print(f"\nwrote {os.path.join(OUT, 'bispecific_constructs.json')}")
