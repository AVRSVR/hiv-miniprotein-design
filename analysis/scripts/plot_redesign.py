"""
Figure 10: the redesign round, plotted as the negative result it is.

Left: dependence on the four variable gp120 positions, stacked per variant, so
the failure is visible as the total bar growing rather than shrinking relative
to the parent. Right: the parent control scored under both protocols, which is
why the per-variant binding energies are not interpretable as affinity changes.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RED = r"C:\Users\anant\Desktop\hiv_miniprotein_design\paper_analysis\07_redesign"
FIG = r"C:\Users\anant\Desktop\hiv_miniprotein_design\paper_analysis\04_figures"
TBL = r"C:\Users\anant\Desktop\hiv_miniprotein_design\paper_analysis\05_tables"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

rows = json.load(open(os.path.join(RED, "redesign_d5.json")))
POS = [("ARG444", "Arg444 (40% cons.)", "#D55E00"),
       ("GLN440", "Gln440 (23% cons.)", "#E69F00"),
       ("ARG308", "Arg308 (28% cons.)", "#CC79A7"),
       ("LYS305", "Lys305 (57% cons.)", "#56B4E9")]
names = [r["variant"] for r in rows]
x = np.arange(len(rows))

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6),
                              gridspec_kw={"width_ratios": [1.9, 1]})

bottom = np.zeros(len(rows))
for key, label, col in POS:
    v = np.array([r["liability_dependence_kcal"].get(key, 0.0) for r in rows])
    ax.bar(x, v, 0.62, bottom=bottom, color=col, label=label,
           edgecolor="white", linewidth=0.6)
    bottom += v
parent = bottom[0]
ax.axhline(parent, ls="--", lw=1.3, color="#444444")
ax.text(len(rows) - 0.45, parent + 0.5, f"parent D5 ({parent:.1f})",
        ha="right", fontsize=9.5, color="#444444")
for xi, tot in zip(x, bottom):
    ax.text(xi, tot + 0.6, f"{tot:.1f}", ha="center", fontsize=9.5)
ax.set_xticks(x)
ax.set_xticklabels([n if n != "WT" else "parent" for n in names])
ax.set_ylabel("Dependence on variable gp120\npositions (kcal/mol)")
ax.set_title("Every mutant increased the liability it was meant to remove",
             fontsize=11.5)
ax.set_ylim(0, bottom.max() * 1.22)
ax.legend(frameon=False, fontsize=9, loc="upper left", ncol=2)
ax.spines[["top", "right"]].set_visible(False)
ax.text(-0.085, 1.0, "a", transform=ax.transAxes, fontsize=17,
        fontweight="bold", va="top", ha="right")

# protocol offset: the same molecule under both protocols
short = rows[0]["dG_kcal_mol"]
short_e = rows[0]["sem_kcal_mol"]
mm = {r["design"]: r for r in json.load(open(os.path.join(TBL, "mmgbsa.json")))}
longs = [mm[t]["dG_bind_kcal_mol"] for t in
         ("Latest_Best", "Latest_Best_rep2", "Latest_Best_rep3") if t in mm]
ax2.bar([0], [short], 0.5, yerr=[short_e], color="#999999",
        edgecolor="black", linewidth=0.7, capsize=5)
ax2.bar([1], [np.mean(longs)], 0.5, yerr=[np.std(longs, ddof=1)],
        color="#009E73", edgecolor="black", linewidth=0.7, capsize=5)
ax2.set_xticks([0, 1])
ax2.set_xticklabels(["0.4 ns\n(redesign round)", "2 ns x3\n(main protocol)"])
ax2.set_ylabel("MM/GBSA $\Delta$G (kcal/mol)")
ax2.set_title("Same molecule, two protocols", fontsize=11.5)
ax2.spines[["top", "right"]].set_visible(False)
ax2.annotate("", xy=(0.5, short), xytext=(0.5, np.mean(longs)),
             arrowprops=dict(arrowstyle="<->", color="#B22222", lw=1.4))
ax2.text(0.56, (short + np.mean(longs)) / 2,
         f"{abs(short - np.mean(longs)):.0f} kcal/mol\noffset",
         fontsize=9.5, color="#B22222", va="center")
ax2.text(-0.22, 1.0, "b", transform=ax2.transAxes, fontsize=17,
         fontweight="bold", va="top", ha="right")

fig.tight_layout()
p = os.path.join(FIG, "fig10_redesign.png")
fig.savefig(p, dpi=300, facecolor="white")
print("wrote", p)
