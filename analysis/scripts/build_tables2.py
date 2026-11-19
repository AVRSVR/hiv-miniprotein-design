"""
Tables 7-9: interface hotspots, escape risk, developability.

Also emits a corrected panel summary. Chai1_Best and model_4_chai carry the
identical 73-residue sequence and are two predicted poses of one molecule, not two
designs; the panel therefore contains four unique sequences.
"""
import os, json, glob
import numpy as np
import pandas as pd

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
TBL = os.path.join(PROJECT, "paper_analysis", "05_tables")

AA3TO1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
          "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
          "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
          "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
NAME = {"S310": "D1 (CD4-mimetic control)", "Chai1_Best": "D2 (de novo bispecific)",
        "Cyclic_Ring": "D3 (CD4-mimetic cyclic)",
        "model_4_chai": "D2' (alternative pose)",
        "Latest_Best": "D5 (de novo, coreceptor site)"}


def one(r):
    return AA3TO1.get(r.upper(), "X")


def save(name, df):
    if df.empty:
        print(f"[{name}] empty")
        return
    df.to_csv(os.path.join(TBL, name + ".csv"), index=False)
    with open(os.path.join(TBL, name + ".md"), "w", encoding="utf-8") as fh:
        fh.write(df.to_markdown(index=False))
    print(f"\n=== {name} ===")
    print(df.to_string(index=False))


# ------------------------------------------------------------ Table 7
rows = []
for p in sorted(glob.glob(os.path.join(TBL, "alascan_*.json"))):
    d = json.load(open(p))
    res = d["residues"]
    hs = [r for r in res if r["hotspot"]]
    top = sorted(hs, key=lambda z: -z["ddG_kcal_mol"])[:5]
    rows.append({
        "Design": NAME.get(d["design"], d["design"]),
        "Interface residues scanned": len(res),
        "Hotspots (>=1 kcal/mol)": len(hs),
        "Hotspot fraction (%)": round(100 * len(hs) / max(len(res), 1), 1),
        "Sum of hotspot ddG (kcal/mol)": round(sum(z["ddG_kcal_mol"] for z in hs), 1),
        "Strongest hotspots": ", ".join(
            f"{one(z['resname'])}{z['resid']} ({z['ddG_kcal_mol']:+.1f})" for z in top),
    })
save("table7_alanine_scan", pd.DataFrame(rows).sort_values(
    "Hotspots (>=1 kcal/mol)", ascending=False))

# ------------------------------------------------------------ Table 8
rows = []
for p in sorted(glob.glob(os.path.join(TBL, "escape_risk_*.json"))):
    d = json.load(open(p))
    cs = [c for c in d["contacts"] if c["conservation"] is not None]
    dep = [c for c in cs if c["ddG_kcal_mol"] >= 1.0]
    resistant = [c for c in dep if c["conservation"] >= 0.95]
    vulnerable = [c for c in dep if c["conservation"] < 0.80]
    top = d.get("top_escape_risk", [])[:3]
    rows.append({
        "Design": NAME.get(d["design"], d["design"]),
        "gp120 contacts scored": len(cs),
        "Load-bearing (>=1 kcal/mol)": len(dep),
        "Escape-resistant (>=95% cons.)": len(resistant),
        "On variable positions (<80% cons.)": len(vulnerable),
        "Escape-resistant fraction (%)": round(100 * len(resistant) / max(len(dep), 1), 1),
        "Highest-risk contacts": ", ".join(
            f"{one(c['resname'])}{c['hxb2_pos']} ({100*c['conservation']:.0f}% cons.)"
            for c in top),
    })
save("table8_escape_risk", pd.DataFrame(rows).sort_values(
    "Escape-resistant fraction (%)", ascending=False))

# ------------------------------------------------------------ Table 9
p = os.path.join(TBL, "developability.json")
if os.path.exists(p):
    d = json.load(open(p))
    seen, rows = set(), []
    for r in d:
        if r["sequence"] in seen:          # D2/D2' duplicate
            continue
        seen.add(r["sequence"])
        rows.append({
            "Design": r["design"], "Length": r["length"],
            "MW (kDa)": r["MW_kDa"], "pI": r["pI"],
            "Net charge pH 7.4": r["net_charge_pH7.4"],
            "GRAVY": r["GRAVY"],
            "Instability index": r["instability_index"],
            "Predicted stable": "yes" if r["stable_by_instability_index"] else "no",
            "Exposed hydrophobic (frac)": r["exposed_hydrophobic_frac"],
            "Largest hydrophobic patch (res)": r["largest_hydrophobic_patch_res"],
            "Cys": r["n_cys"], "Unpaired Cys": r["n_free_cys_if_unpaired"],
            "N-glyc sequons": r["n_glycosylation_sequons"],
            "Deamidation motifs": r["n_deamidation_NG_NS"],
            "Helix/sheet (%)": f"{100*r['helix_frac']:.0f}/{100*r['sheet_frac']:.0f}",
        })
    save("table9_developability", pd.DataFrame(rows))
