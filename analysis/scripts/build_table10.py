"""Table 10: D5 redesign round against the escape liability."""
import os, json
import pandas as pd

TBL = r"C:\Users\anant\Desktop\hiv_miniprotein_design\paper_analysis\05_tables"
RED = r"C:\Users\anant\Desktop\hiv_miniprotein_design\paper_analysis\07_redesign"

rows = []
for r in json.load(open(os.path.join(RED, "redesign_d5.json"))):
    dep = r["liability_dependence_kcal"]
    rows.append({
        "Variant": r["variant"],
        "Mutation(s)": ", ".join(r["mutations"]) or "none (parent D5)",
        "MM/GBSA dG (kcal/mol)": f"{r['dG_kcal_mol']:.1f} +/- {r['sem_kcal_mol']:.1f}",
        "Arg444 (36% cons.)": dep.get("ARG444"),
        "Gln440 (23% cons.)": dep.get("GLN440"),
        "Arg308 (37% cons.)": dep.get("ARG308"),
        "Lys305 (68% cons.)": dep.get("LYS305"),
        "Total liability (kcal/mol)": r["liability_total_kcal"],
        "Meets criterion": "no",      # binding held but liability rose in every case
    })
df = pd.DataFrame(rows)
df.to_csv(os.path.join(TBL, "table10_redesign.csv"), index=False)
with open(os.path.join(TBL, "table10_redesign.md"), "w", encoding="utf-8") as fh:
    fh.write(df.to_markdown(index=False))
print(df.to_string(index=False))

# ------------------------------------------------- Table 12: footprint overlap
fp = os.path.join(TBL, "footprint_overlap.json")
if os.path.exists(fp):
    d = json.load(open(fp))
    df12 = pd.DataFrame([{
        "Design": r["label"],
        "Binder footprint (res)": r["binder_footprint_res"],
        "Shared with CD4 footprint": f"{r['shared_with_CD4']}/{r['cd4_footprint_res']}",
        "% of CD4 footprint": r["pct_of_CD4_footprint"],
        "Shared with 17b footprint": f"{r['shared_with_17b']}/{r['17b_footprint_res']}",
        "% of 17b footprint": r["pct_of_17b_footprint"],
    } for r in d["per_design"]])
    df12.to_csv(os.path.join(TBL, "table12_footprint_overlap.csv"), index=False)
    with open(os.path.join(TBL, "table12_footprint_overlap.md"), "w",
              encoding="utf-8") as fh:
        fh.write(df12.to_markdown(index=False))
    print("\n=== table12_footprint_overlap ===")
    print(df12.to_string(index=False))
