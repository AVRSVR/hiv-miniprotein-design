"""
RFdiffusion vs Rosetta design-arm comparison with significance testing.

AlphaFold2 confidence metrics (pLDDT, mean PAE, pTM) for both arms, compared with
Mann-Whitney U (non-parametric; these distributions are skewed and unequal-n) plus
rank-biserial correlation as effect size. Also summarizes the off-target (BLAST vs
human proteome) and MHC-I immunogenicity screens.
"""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")
os.makedirs(OUT, exist_ok=True)

ros = pd.read_csv(os.path.join(PROJECT, "af2_rosetta_scores.csv"))
rfd = pd.read_csv(os.path.join(PROJECT, "af2_rfdiffusion_all_scores.csv"))

ros = ros.rename(columns={"mean_pae": "pae"})
for df in (ros, rfd):
    for c in ("plddt", "pae", "ptm"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

ros_c = ros.dropna(subset=["plddt", "pae", "ptm"])
rfd_c = rfd.dropna(subset=["plddt", "pae", "ptm"])


def rank_biserial(a, b):
    """Effect size from Mann-Whitney U: 2*U/(n1*n2) - 1, range [-1,1]."""
    u = mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2.0 * u / (len(a) * len(b)) - 1.0


def desc(x):
    x = np.asarray(x, float)
    return {"n": int(x.size), "mean": round(float(x.mean()), 2),
            "sd": round(float(x.std(ddof=1)), 2), "median": round(float(np.median(x)), 2),
            "min": round(float(x.min()), 2), "max": round(float(x.max()), 2)}


rows = []
for metric, better in (("plddt", "higher"), ("pae", "lower"), ("ptm", "higher")):
    a, b = ros_c[metric].values, rfd_c[metric].values
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    rows.append({
        "metric": metric, "better": better,
        "rosetta": desc(a), "rfdiffusion": desc(b),
        "mannwhitney_U": float(u), "p_value": float(p),
        "rank_biserial": round(rank_biserial(a, b), 3),
        "significant_p<0.001": bool(p < 1e-3),
    })

# AF2 confidence pass rates (pLDDT>=80 & PAE<=5 is a common binder-design bar)
def pass_rate(df):
    ok = (df["plddt"] >= 80) & (df["pae"] <= 5)
    return {"n": int(len(df)), "n_pass": int(ok.sum()),
            "pct": round(100.0 * ok.mean(), 1)}

summary = {
    "af2_metrics": rows,
    "high_confidence_rate_plddt80_pae5": {
        "rosetta": pass_rate(ros_c), "rfdiffusion": pass_rate(rfd_c),
    },
}

# screening summaries
bl_p = os.path.join(PROJECT, "blast_screening_results.csv")
if os.path.exists(bl_p):
    bl = pd.read_csv(bl_p)
    col = "flag" if "flag" in bl.columns else None
    arm = "arm" if "arm" in bl.columns else None
    if col and arm:
        summary["blast_human_offtarget"] = {
            str(k): {str(kk): int(vv) for kk, vv in v.items()}
            for k, v in bl.groupby(arm)[col].value_counts().unstack(fill_value=0).iterrows()
        }

mh_p = os.path.join(PROJECT, "mhcflurry_results.csv")
if os.path.exists(mh_p):
    mh = pd.read_csv(mh_p)
    fcol = [c for c in mh.columns if "flag" in c.lower()]
    icol = [c for c in mh.columns if "ic50" in c.lower()]
    d = {"n": int(len(mh))}
    if fcol:
        d["flag_counts"] = {str(k): int(v) for k, v in mh[fcol[0]].value_counts().items()}
    if icol:
        v = pd.to_numeric(mh[icol[0]], errors="coerce").dropna()
        if len(v):
            d["best_ic50_nM"] = {"median": round(float(v.median()), 1),
                                 "min": round(float(v.min()), 1),
                                 "n_strong_binders_lt50nM": int((v < 50).sum())}
    summary["mhc_i_immunogenicity"] = d

with open(os.path.join(OUT, "arm_comparison.json"), "w") as f:
    json.dump(summary, f, indent=2)

print("=== AlphaFold2 confidence: Rosetta vs RFdiffusion ===")
print(f"{'metric':<8}{'Rosetta mean+-SD':>22}{'RFdiff mean+-SD':>22}{'p':>12}{'effect':>9}")
for r in rows:
    R, F = r["rosetta"], r["rfdiffusion"]
    print(f"{r['metric']:<8}{R['mean']:>13.2f} ± {R['sd']:<6.2f}{F['mean']:>13.2f} ± {F['sd']:<6.2f}"
          f"{r['p_value']:>12.3e}{r['rank_biserial']:>9.3f}")
print(f"\nn: Rosetta={len(ros_c)}, RFdiffusion={len(rfd_c)}")
hc = summary["high_confidence_rate_plddt80_pae5"]
print(f"High-confidence (pLDDT>=80 & PAE<=5): "
      f"Rosetta {hc['rosetta']['n_pass']}/{hc['rosetta']['n']} ({hc['rosetta']['pct']}%)  |  "
      f"RFdiffusion {hc['rfdiffusion']['n_pass']}/{hc['rfdiffusion']['n']} ({hc['rfdiffusion']['pct']}%)")
for k in ("blast_human_offtarget", "mhc_i_immunogenicity"):
    if k in summary:
        print(f"\n{k}: {json.dumps(summary[k])}")
