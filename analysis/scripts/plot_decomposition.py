"""
Figure 11: is the designed interface electrostatics-driven?

Each bar splits the MM/GBSA binding energy into its van der Waals part and
everything electrostatic (Coulomb plus polar solvation). The experimental
CD4-gp120 and b12-gp120 interfaces are scored identically and bracket the range
a real gp120 interface occupies.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PA = r"C:\Users\anant\Desktop\hiv_miniprotein_design\paper_analysis"
rows = json.load(open(os.path.join(PA, "05_tables", "energy_decomposition.json")))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

order = ["D1 (CD4-mimetic control)", "D2 (de novo bispecific)",
         "D3 (CD4-mimetic cyclic)", "D2' (alternative pose)",
         "D5 (de novo, coreceptor site)", "CD4-gp120 (experimental)",
         "b12-gp120 (experimental)"]
rows = sorted(rows, key=lambda r: order.index(r["system"]) if r["system"] in order else 99)
lab = [r["system"].split(" (")[0] + ("\n(experimental)" if "experimental" in r["system"] else "")
       for r in rows]
vdw = np.array([r["dG_vdW_kcal"] for r in rows])
ele = np.array([r["dG_elec_solv_kcal"] for r in rows])
tot = np.array([r["dG_total_kcal"] for r in rows])
exp = np.array(["experimental" in r["system"] for r in rows])
x = np.arange(len(rows))

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.8),
                              gridspec_kw={"width_ratios": [1.65, 1]})
ax.bar(x, vdw, 0.6, color="#0072B2", label="van der Waals")
ax.bar(x, ele, 0.6, color="#E69F00", label="Electrostatic + polar solvation")
ax.plot(x, tot, "o", color="black", ms=7, label="Net binding energy", zorder=5)
ax.axhline(0, color="black", lw=0.9)
for xi, t in zip(x, tot):
    ax.annotate(f"{t:.0f}", (xi, t), textcoords="offset points", xytext=(0, 9),
                ha="center", fontsize=9)
ax.set_xticks(x)
ax.set_xticklabels(lab, fontsize=9)
ax.set_ylabel("Energy component (kcal/mol)")
ax.set_title("Binding is van der Waals-driven in every interface, designed or not",
             fontsize=11.5)
ax.legend(frameon=False, fontsize=9, loc="lower left")
ax.spines[["top", "right"]].set_visible(False)
ax.text(-0.075, 1.0, "a", transform=ax.transAxes, fontsize=17,
        fontweight="bold", va="top", ha="right")

frac = -ele / -tot * -1          # elec contribution as a fraction of the net
frac = np.array([r["elec_fraction"] for r in rows])
cols = ["#999999" if e else "#009E73" for e in exp]
cols = [("#D55E00" if "D5" in r["system"] else ("#555555" if e else "#9ECAE1"))
        for r, e in zip(rows, exp)]
ax2.barh(x, frac, 0.6, color=cols)
lo, hi = frac[exp].min(), frac[exp].max()
ax2.axvspan(lo, hi, color="#000000", alpha=0.07, zorder=0)
ax2.text((lo + hi) / 2, len(rows) - 0.3, "range spanned by the\nexperimental complexes",
         ha="center", fontsize=8.5, color="#444444")
ax2.set_yticks(x)
ax2.set_yticklabels(lab, fontsize=9)
ax2.invert_yaxis()
ax2.set_xlabel("Electrostatic fraction of net binding energy")
ax2.set_title("D5 falls inside that range", fontsize=11.5)
ax2.spines[["top", "right"]].set_visible(False)
ax2.text(-0.30, 1.0, "b", transform=ax2.transAxes, fontsize=17,
         fontweight="bold", va="top", ha="right")

fig.tight_layout()
p = os.path.join(PA, "04_figures", "fig11_energy_decomposition.png")
fig.savefig(p, dpi=300, facecolor="white")
print("wrote", p)

# Table 13
import pandas as pd
df = pd.DataFrame([{
    "Interface": r["system"],
    "Net dG (kcal/mol)": r["dG_total_kcal"],
    "van der Waals": r["dG_vdW_kcal"],
    "Electrostatic + polar solvation": r["dG_elec_solv_kcal"],
    "Electrostatic fraction": r["elec_fraction"],
    "Frames": r["n_frames"],
} for r in rows])
T = os.path.join(PA, "05_tables")
df.to_csv(os.path.join(T, "table13_energy_decomposition.csv"), index=False)
open(os.path.join(T, "table13_energy_decomposition.md"), "w",
     encoding="utf-8").write(df.to_markdown(index=False))
print(df.to_string(index=False))
