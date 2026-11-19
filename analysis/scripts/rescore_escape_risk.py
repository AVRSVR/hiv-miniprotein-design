"""
Re-attach conservation to the existing escape-risk scans.

The escape score has two inputs: a per-contact Delta-Delta-G from the decoupling scan,
and the conservation of that gp120 position. Only the second depends on which Env
sequences were used, so expanding the sequence set does not require re-running the
scan - which would mean re-opening every trajectory and rebuilding OpenMM contexts
for a result that cannot change. The stored Delta-Delta-G values are reused and
conservation, the escape score and the derived summaries are recomputed from the
shared cache.

Run: python rescore_escape_risk.py
"""
import glob
import io
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_core

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "05_tables")
RESISTANT = 0.95
LOAD_BEARING = 1.0

cons, meta = conservation_core.load()
print("conservation: %d sequences, %d clusters at %.0f%% identity, %d positions"
      % (meta["sequences_used"], meta["clusters_at_98pct_identity"],
         100 * meta["cluster_identity"], meta["positions_scored"]))

for path in sorted(glob.glob(os.path.join(OUT, "escape_risk_*.json"))):
    with io.open(path, encoding="utf-8") as fh:
        d = json.load(fh)

    for c in d["contacts"]:
        v = cons.get(c["hxb2_pos"])
        c["conservation"] = round(v, 3) if v is not None else None

    scored = [c for c in d["contacts"]
              if c["conservation"] is not None and c["ddG_kcal_mol"] > 0]
    mx = max((c["ddG_kcal_mol"] for c in scored), default=0.0)
    for c in d["contacts"]:
        if c["conservation"] is None or c["ddG_kcal_mol"] <= 0 or mx <= 0:
            c["escape_risk"] = None
            continue
        c["escape_risk"] = round((c["ddG_kcal_mol"] / mx)
                                 * (1.0 - c["conservation"]), 3)

    load_bearing = [c for c in d["contacts"]
                    if c["conservation"] is not None
                    and c["ddG_kcal_mol"] >= LOAD_BEARING]
    resistant = [c for c in load_bearing if c["conservation"] >= RESISTANT]
    d["top_escape_risk"] = sorted(
        [c for c in d["contacts"] if c.get("escape_risk") is not None],
        key=lambda z: -z["escape_risk"])[:10]
    d["escape_resistant_contacts"] = sorted(resistant,
                                            key=lambda z: -z["ddG_kcal_mol"])
    d["conservation_source"] = {
        "sequences": meta["sequences_used"],
        "clusters": meta["clusters_at_98pct_identity"],
        "estimator": "modal residue frequency over 98%-identity cluster representatives",
    }

    with io.open(path, "w", encoding="utf-8") as fh:
        json.dump(d, fh, indent=2)

    pct = 100.0 * len(resistant) / len(load_bearing) if load_bearing else 0.0
    below80 = sum(1 for c in load_bearing if c["conservation"] < 0.80)
    print("\n%-4s load-bearing %2d | escape-resistant %2d (%.1f%%) | below 80%%: %d"
          % (d["label"], len(load_bearing), len(resistant), pct, below80))
    for c in sorted(load_bearing, key=lambda z: -z["ddG_kcal_mol"])[:5]:
        print("   %s%-4s ddG %+6.2f  cons %3.0f%%  risk %s"
              % (c["resname"], c["hxb2_pos"], c["ddG_kcal_mol"],
                 100 * c["conservation"], c["escape_risk"]))
print("\nrewrote %d scans" % len(glob.glob(os.path.join(OUT, "escape_risk_*.json"))))
