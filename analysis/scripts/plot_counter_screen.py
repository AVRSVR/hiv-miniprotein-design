"""
Figure 13: the counter-screen.

One point per minimised face-on placement of D5, four targets screened identically, with
the designed pose marked. The point of the figure is the gap between the designed pose
and everything else, and the fact that arbitrary contact scores about the same on every
surface including the target's.
"""
import io
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TBL = os.path.join(PA, "05_tables")
FIG = os.path.join(PA, "04_figures")

INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#dcdcdc"
NAME = {"gp120": "gp120\n(the target)", "HSA_1AO6": "serum albumin\n1AO6",
        "lysozyme_1DPX": "lysozyme\n1DPX", "ubiquitin_1UBQ": "ubiquitin\n1UBQ"}
COLOR = {"gp120": "#009E73", "HSA_1AO6": "#E69F00",
         "lysozyme_1DPX": "#CC79A7", "ubiquitin_1UBQ": "#56B4E9"}

plt.rcParams.update({"font.family": "Arial", "font.size": 9,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.linewidth": 0.8, "figure.dpi": 300})

d = json.load(io.open(os.path.join(TBL, "counter_screen.json"), encoding="utf-8"))
dz = d["designs"]["Latest_Best"]
designed = dz["designed_pose_on_gp120_kcal"]
order = [k for k in ("gp120", "HSA_1AO6", "lysozyme_1DPX", "ubiquitin_1UBQ")
         if k in dz["targets"]]

rng = np.random.default_rng(7)
fig, ax = plt.subplots(figsize=(6.6, 4.4))
for i, k in enumerate(order):
    t = dz["targets"][k]
    # the stored summary carries best, median, mean and sd; draw the span and the marks
    ax.plot([i, i], [t["best_kcal"], t["best_kcal"] + 2 * (t["sd_kcal"] or 0)],
            color=COLOR[k], lw=1.0, alpha=0.5, zorder=1)
    ax.scatter([i], [t["best_kcal"]], s=70, color=COLOR[k], zorder=3,
               edgecolor="white", linewidth=0.9, label=None)
    ax.scatter([i], [t["median_kcal"]], s=70, marker="D", facecolor="white",
               edgecolor=COLOR[k], linewidth=1.6, zorder=3)
    ax.text(i, t["best_kcal"] - 1.6, "%.1f" % t["best_kcal"], ha="center",
            fontsize=7.5, color=INK)

ax.axhline(designed, color="#009E73", lw=1.2, ls=(0, (4, 3)), zorder=2)
ax.text(len(order) - 0.52, designed + 1.4,
        "D5 designed pose on gp120, %.1f kcal/mol" % designed,
        ha="right", fontsize=8, color="#097a5c")

ax.set_xticks(range(len(order)))
ax.set_xticklabels([NAME[k] for k in order], fontsize=8.5)
ax.set_ylabel("MM/GBSA of the placement (kcal/mol)")
ax.set_title("Face-on placements scored identically on four surfaces",
             loc="left", fontsize=10, color=INK, pad=8)
ax.grid(True, axis="y", color=GRID, lw=0.6)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.set_xlim(-0.6, len(order) - 0.4)
ax.set_ylim(designed - 6, 2)

from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", color="none", markerfacecolor=MUTED,
                          markersize=8, label="best of 8 placements"),
                   Line2D([], [], marker="D", color="none", markeredgecolor=MUTED,
                          markerfacecolor="white", markersize=8, label="median")],
          frameon=False, fontsize=8, loc="center left", labelcolor=INK)
fig.tight_layout()
out = os.path.join(FIG, "fig13_counter_screen.png")
fig.savefig(out, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote", os.path.basename(out))
