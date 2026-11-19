"""
Assemble publication tables from the completed analyses.

Table 1  design panel: size, sequence, interface metrics, epitope class
Table 2  CD4 competition and 17b overlap
Table 3  pairwise co-binding compatibility on gp120
Table 4  MD interface stability (whatever runs have finished)
"""
import os, json, glob
import numpy as np
import pandas as pd

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SDIR = os.path.join(PROJECT, "paper_analysis", "01_structures")
MDIR = os.path.join(PROJECT, "paper_analysis", "03_md")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")
os.makedirs(OUT, exist_ok=True)

# panel labels, matching the manuscript
LABEL = {
    "S310": "D1 (CD4-mimetic control)",
    "Chai1_Best": "D2 (de novo bispecific)",
    "Cyclic_Ring": "D3 (CD4-mimetic cyclic)",
    "model_4_chai": "D2' (alternative pose)",
    "Latest_Best": "D5 (de novo, coreceptor site)",
}
# origin, from pairwise alignment against human CD4 (UniProt P01730)
PROVENANCE = {
    "S310": "CD4-mimetic (52.6% id to CD4 D1 over 97 aa)",
    "Cyclic_Ring": "CD4-mimetic (100% id over 11 aa)",
    "Chai1_Best": "de novo",
    "model_4_chai": "de novo (same sequence as D2)",
    "Latest_Best": "de novo",
}
EXCLUDE = {"docked_complex"}   # superposition failed (>10 A); not interpretable

iface = {r["design"]: r for r in
         json.load(open(os.path.join(SDIR, "interface_metrics.json")))}
epi = {r["design"]: r for r in
       json.load(open(os.path.join(SDIR, "epitope_verification.json")))}
comp_all = json.load(open(os.path.join(OUT, "cd4_competition.json")))
comp = {r["design"]: r for r in comp_all["per_design"]}

designs = [d for d in iface if d not in EXCLUDE]
order = [d for d in LABEL if d in designs] + [d for d in designs if d not in LABEL]

# ---------------------------------------------------------------- Table 1
t1 = []
for d in order:
    i, e = iface[d], epi.get(d, {})
    t1.append({
        "Design": LABEL.get(d, d),
        "Origin": PROVENANCE.get(d, "unknown"),
        "Length (aa)": i["binder_len"],
        "Contacts (<4.5 A)": i["n_contacts"],
        "H-bonds": i["n_hbonds"],
        "Salt bridges": i["n_saltbridges"],
        "BSA (A^2)": i["dSASA_A2"],
        "Interface res (gp120)": i["n_iface_res_receptor"],
        "Interface res (binder)": i["n_iface_res_binder"],
        "CD4bs hotspots engaged": f"{i['n_hotspots_engaged']}/14",
        "Epitope class": ("CD4-binding site" if e.get("verdict") == "CD4bs-engaged"
                          else "CD4bs-peripheral" if e.get("verdict") == "peripheral"
                          else "Non-CD4bs"),
    })
df1 = pd.DataFrame(t1)

# ---------------------------------------------------------------- Table 2
t2 = []
for d in order:
    c = comp.get(d)
    if not c:
        continue
    t2.append({
        "Design": LABEL.get(d, d),
        "Fit RMSD to 1GC1 (A)": c["cd4_epitope_fit_RMSD_A"],
        "Min dist to CD4 (A)": c["binder_CD4_min_dist_A"],
        "CD4 clash atoms (<3 A)": c["binder_CD4_clash_atoms_3A"],
        "CD4 overlap atoms (<4.5 A)": c["binder_CD4_contact_atoms_4.5A"],
        "Min dist to 17b (A)": c["binder_17b_min_dist_A"],
        "CD4 competition": c["CD4_competition"],
    })
df2 = pd.DataFrame(t2)

# ---------------------------------------------------------------- Table 3
t3 = []
for p in comp_all["pairwise_cobinding"]:
    if p["design_1"] in EXCLUDE or p["design_2"] in EXCLUDE:
        continue
    t3.append({
        "Design A": LABEL.get(p["design_1"], p["design_1"]),
        "Design B": LABEL.get(p["design_2"], p["design_2"]),
        "Min dist (A)": p["min_dist_A"],
        "Clashing atoms (<3 A)": p["clash_atoms_3A"],
        "Co-binding": p["co_binding"],
    })
df3 = pd.DataFrame(t3).sort_values("Clashing atoms (<3 A)")

# ---------------------------------------------------------------- Table 4
t4 = []
for f in sorted(glob.glob(os.path.join(MDIR, "*_md_metrics.json"))):
    m = json.load(open(f))
    d = m["design"]
    if d in EXCLUDE or d.startswith(("TEST", "SPEED")):
        continue
    fr = m["frames"]
    if not fr:
        continue
    ir = np.array([x["iRMSD_A"] for x in fr])
    br = np.array([x["binder_RMSD_A"] for x in fr])
    ct = np.array([x["contacts"] for x in fr])
    md = np.array([x["min_dist_A"] for x in fr])
    half = len(ir) // 2
    t4.append({
        "_tag": d,
        "Design": LABEL.get(d, d),
        "Sim length (ns)": m["ns"],
        "iRMSD mean (A)": round(float(ir.mean()), 2),
        "iRMSD final-half mean (A)": round(float(ir[half:].mean()), 2),
        "iRMSD max (A)": round(float(ir.max()), 2),
        "Binder fold RMSD (A)": round(float(br.mean()), 2),
        "Contacts mean": int(ct.mean()),
        # baseline is the mean of the opening frames, not a single noisy frame
        "Contacts retained (%)": round(
            100.0 * ct[half:].mean() / max(ct[:max(3, len(ct) // 10)].mean(), 1), 1),
        "Min dist mean (A)": round(float(md.mean()), 2),
        "Bound at end": "yes" if md[-1] < 5.0 else "no",
    })
df4 = pd.DataFrame(t4)

# MM/GBSA binding energies, for whichever trajectories have been processed
mm_p = os.path.join(OUT, "mmgbsa.json")
if os.path.exists(mm_p) and not df4.empty:
    mm = {r["design"]: r for r in json.load(open(mm_p))}
    df4["MM/GBSA dG (kcal/mol)"] = df4["_tag"].map(
        lambda t: mm.get(t, {}).get("dG_bind_kcal_mol"))
    df4["MM/GBSA SEM"] = df4["_tag"].map(
        lambda t: mm.get(t, {}).get("sem_kcal_mol"))
if "_tag" in df4.columns:
    df4 = df4.drop(columns=["_tag"])

# ------------------------------------------------- Table 5: interface benchmark
bench_p = os.path.join(OUT, "interface_benchmark.json")
df5 = pd.DataFrame()
if os.path.exists(bench_p):
    b = json.load(open(bench_p))
    df5 = pd.DataFrame([{
        "Interface": r["interface"], "Type": r["type"],
        "Contacts (<4.5 A)": r["n_contacts"], "H-bonds": r["n_hbonds"],
        "Salt bridges": r["n_saltbridges"], "BSA (A^2)": r["dSASA_A2"],
        "Interface res (gp120)": r["n_iface_res_gp120"],
        "Interface res (partner)": r["n_iface_res_partner"],
    } for r in b]).sort_values(["Type", "BSA (A^2)"], ascending=[True, False])

# ------------------------------------------------- Table 6: bispecific constructs
bis_p = os.path.join(OUT, "bispecific_constructs.json")
df6 = pd.DataFrame()
if os.path.exists(bis_p):
    bb = json.load(open(bis_p))["constructs"]
    df6 = pd.DataFrame([{
        "Construct (N->C)": r["fusion"], "Arm 1": r["arm_1"], "Arm 2": r["arm_2"],
        "C-to-N distance (A)": r["C_to_N_distance_A"],
        "Min linker (res)": r["min_linker_residues"],
        "Suggested linker (res)": r["suggested_linker_residues"],
        "Feasible": "yes" if r["feasible"] else "no",
    } for r in bb])

# ------------------------------------------------- S1: epitope conservation
cons_p = os.path.join(OUT, "epitope_conservation.json")
dfS1 = pd.DataFrame()
if os.path.exists(cons_p):
    c = json.load(open(cons_p))
    dfS1 = pd.DataFrame(c["per_design"])
    dfS1.attrs["gp120_mean"] = c["gp120_mean_conservation_pct"]

for name, df in (("table1_design_panel", df1), ("table2_cd4_competition", df2),
                 ("table3_cobinding", df3), ("table4_md_stability", df4),
                 ("table5_interface_benchmark", df5),
                 ("table6_bispecific_constructs", df6),
                 ("tableS1_epitope_conservation", dfS1)):
    if df.empty:
        print(f"[{name}] empty - skipped")
        continue
    df.to_csv(os.path.join(OUT, name + ".csv"), index=False)
    with open(os.path.join(OUT, name + ".md"), "w") as fh:
        fh.write(df.to_markdown(index=False))
    print(f"\n=== {name} ===")
    print(df.to_string(index=False))
