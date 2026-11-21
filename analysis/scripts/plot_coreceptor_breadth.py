"""
Figure 14: the coreceptor itself, and how much of the epidemic presents each interface.

(a) clashes with CCR5 carried into the open trimer, split N-terminal segment / body
(b) share of CCR5's gp120 footprint each partner covers
(c) median interface retained across 2,948 representatives, against a size-matched null
(d) for every D5 contact, the frequency of the residue the design holds beside the modal
    one, which is the distinction the whole section rests on

Run: python plot_coreceptor_breadth.py
"""
import io
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TBL = os.path.join(PA, "05_tables")
FIG = os.path.join(PA, "04_figures")

INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#dcdcdc"
BLUE, ORANGE, GREEN = "#0072B2", "#E69F00", "#009E73"
VERM, PURPLE, GREY = "#D55E00", "#CC79A7", "#9a9a9a"
COL = {"CD4": "#8bb8d8", "17b": ORANGE, "D5": GREEN, "D2": ORANGE,
       "D2p": VERM, "D3": PURPLE}

plt.rcParams.update({"font.family": "Arial", "font.size": 9,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.linewidth": 0.8, "figure.dpi": 300})


def load(n):
    with io.open(os.path.join(TBL, n), encoding="utf-8") as fh:
        return json.load(fh)


cc = load("ccr5_competition.json")
sc = load("strain_coverage.json")

fig, ax = plt.subplots(2, 2, figsize=(10.4, 7.6))
for a in ax.ravel():
    a.grid(True, axis="y", color=GRID, lw=0.6)
    a.set_axisbelow(True)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)

ORDER = ["CD4", "17b", "D5", "D2", "D2p", "D3"]
NAME = {"CD4": "CD4", "17b": "17b", "D5": "D5", "D2": "D2", "D2p": "D2′", "D3": "D3"}

# ---------------------------------------------------------------- (a) clashes with CCR5
A = ax[0, 0]
x = np.arange(len(ORDER))
nt = [cc["partners"][k]["nterm_clashes_3A"] for k in ORDER]
bd = [cc["partners"][k]["body_clashes_3A"] for k in ORDER]
A.bar(x - 0.19, np.maximum(nt, 0.4), 0.36, color=[COL[k] for k in ORDER],
      edgecolor="white", lw=0.6)
A.bar(x + 0.19, np.maximum(bd, 0.4), 0.36, color=[COL[k] for k in ORDER],
      edgecolor="white", lw=0.6, hatch="////", alpha=0.75)
A.set_yscale("log")
A.set_ylim(0.3, 900)
for xi, (n, b) in enumerate(zip(nt, bd)):
    A.text(xi - 0.19, max(n, 0.4) * 1.25, str(n), ha="center", fontsize=7, color=INK)
    A.text(xi + 0.19, max(b, 0.4) * 1.25, str(b), ha="center", fontsize=7, color=INK)
A.set_xticks(x)
A.set_xticklabels([NAME[k] for k in ORDER])
A.set_ylabel("atom pairs below 3 Å (log)")
A.set_title("a   Steric conflict with CCR5 in the open trimer", loc="left",
            fontsize=10, color=INK, pad=8)
A.legend(handles=[Patch(facecolor=GREY, edgecolor="white",
                        label="N-terminal segment (1-30)"),
                  Patch(facecolor=GREY, edgecolor="white", hatch="////", alpha=0.75,
                        label="receptor body")],
         frameon=False, fontsize=8, loc="upper right", labelcolor=INK)
A.annotate("must not clash", xy=(0, 0.6), xytext=(0.35, 6), fontsize=7.5, color="#2b6ca3",
           arrowprops=dict(arrowstyle="-", lw=0.7, color="#2b6ca3"))

# ------------------------------------------------------- (b) share of CCR5's footprint
B = ax[0, 1]
cov = [cc["partners"][k]["pct_of_ccr5_footprint_covered"] for k in ORDER]
B.bar(x, cov, 0.6, color=[COL[k] for k in ORDER], edgecolor="white", lw=0.8)
for xi, v in zip(x, cov):
    B.text(xi, v + 1.1, "%.1f%%" % v, ha="center", fontsize=7.5, color=INK)
B.set_xticks(x)
B.set_xticklabels([NAME[k] for k in ORDER])
B.set_ylim(0, 42)
B.set_ylabel("share of the CCR5 footprint covered (%)")
B.set_title("b   How much of the receptor's own gp120 footprint", loc="left",
            fontsize=10, color=INK, pad=8)
B.axhline(cc["partners"]["17b"]["pct_of_ccr5_footprint_covered"], color=ORANGE,
          lw=1.0, ls=(0, (4, 3)), zorder=1)
B.text(len(ORDER) - 0.45, cc["partners"]["17b"]["pct_of_ccr5_footprint_covered"] + 1.2,
       "17b, the usual surrogate", ha="right", fontsize=7.5, color="#a9741a")

# ------------------------------------------------------------- (c) cross-strain coverage
C = ax[1, 0]
SETS = [("D1", "D1"), ("D5", "D5"), ("D2", "D2"),
        ("CD4 epitope", "CD4"), ("17b epitope", "17b"), ("CCR5", "CCR5")]
SCOL = {"D1": BLUE, "D5": GREEN, "D2": ORANGE, "CD4 epitope": "#8bb8d8",
        "17b epitope": "#d8a84a", "CCR5": "#4a7f6e"}
xs = np.arange(len(SETS))
med = [sc["sets"][k]["median_retained_pct"] for k, _ in SETS]
n5 = sc["null_model"]["D5"]
C.axhspan(n5["p5_pct"], n5["p95_pct"], color="#efefef", zorder=0)
C.axhline(n5["median_of_medians_pct"], color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
C.text(-0.42, n5["p95_pct"] + 1.4,
       "size-matched random surface patches: median %.1f%%, 5th to 95th percentile "
       "shaded" % n5["median_of_medians_pct"], ha="left", fontsize=7.5, color=MUTED)
C.bar(xs, med, 0.6, color=[SCOL[k] for k, _ in SETS], edgecolor="white", lw=0.8,
      zorder=2)
for xi, v in zip(xs, med):
    C.text(xi, v + 1.2, "%.1f" % v, ha="center", fontsize=7.5, color=INK, zorder=3)
C.set_xticks(xs)
C.set_xticklabels([lab for _, lab in SETS])
C.set_ylim(0, 100)
C.set_ylabel("median of the interface retained per virus (%)")
C.set_title("c   Coverage across 2,948 cluster representatives", loc="left",
            fontsize=10, color=INK, pad=8)

# ----------------------------------------------- (d) whose residue, position by position
D = ax[1, 1]
pp = sorted(sc["sets"]["D5"]["per_position"], key=lambda q: -q["weight_kcal"])
xp = np.arange(len(pp))
tgt = [100 * q["bg505_frequency"] for q in pp]
mod = [100 * q["modal_frequency"] for q in pp]
D.vlines(xp, tgt, mod, color=GRID, lw=2.2, zorder=1)
D.scatter(xp, mod, s=28, facecolor="white", edgecolor=MUTED, lw=1.1, zorder=3)
D.scatter(xp, tgt, s=30, color=[GREEN if q["bg505_is_modal"] else VERM for q in pp],
          zorder=4)
D.set_xticks(xp)
D.set_xticklabels(["%s%d" % (q["bg505_residue"], q["position"]) for q in pp],
                  rotation=90, fontsize=7)
D.set_ylim(0, 105)
D.set_ylabel("frequency across representatives (%)")
D.set_title("d   D5's contacts: the residue it holds, and the modal one", loc="left",
            fontsize=10, color=INK, pad=8)
D.legend(handles=[Line2D([], [], marker="o", color="none", markerfacecolor=GREEN,
                         markersize=7, label="design holds the modal residue"),
                  Line2D([], [], marker="o", color="none", markerfacecolor=VERM,
                         markersize=7, label="design holds a minority residue"),
                  Line2D([], [], marker="o", color="none", markeredgecolor=MUTED,
                         markerfacecolor="white", markersize=7,
                         label="modal residue at that position")],
         frameon=False, fontsize=7.5, loc="lower left", labelcolor=INK)
D.text(len(pp) - 0.3, 101, "contacts ordered by ΔΔG", ha="right", fontsize=7,
       color=MUTED)

fig.tight_layout(h_pad=2.4, w_pad=2.0)
out = os.path.join(FIG, "fig14_coreceptor_and_breadth.png")
fig.savefig(out, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote", os.path.basename(out))
