"""
Figure 12: the designs measured in native Env context.

(a) occlusion on the closed glycosylated trimer, glycan and neighbouring-subunit
    clashes side by side, with CD4 and 17b as the experimental references
(b) how much of each epitope's surface survives trimer assembly and glycosylation
(c) the share of each design's alanine-scan hotspot energy that sits on residues in
    glycan contact, which is what decides whether the clash matters
(d) whether the footprint keeps its shape through CD4-induced opening, same isolate
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
COLOR = {"D5": "#009E73", "D1": "#0072B2", "D2": "#E69F00",
         "D2p": "#D55E00", "D3": "#CC79A7", "CD4": "#56B4E9", "17b": "#999999"}
SHORT = {"CD4": "CD4", "17b": "17b", "D1": "D1", "D2": "D2", "D2p": "D2'",
         "D3": "D3", "D5": "D5"}

plt.rcParams.update({"font.family": "Arial", "font.size": 9,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.linewidth": 0.8, "figure.dpi": 300})

nc = json.load(io.open(os.path.join(TBL, "native_context.json"), encoding="utf-8"))
sd = json.load(io.open(os.path.join(TBL, "state_dependence.json"), encoding="utf-8"))
P = nc["partners"]

order = ["CD4", "17b", "D1", "D2", "D2p", "D3", "D5"]
fig, axes = plt.subplots(2, 2, figsize=(7.6, 6.4))


def dress(ax, panel, title):
    ax.text(-0.22, 1.16, panel, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top", ha="left")
    ax.set_title(title, loc="left", fontsize=9, color=INK, pad=6)
    ax.grid(True, axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# --------------------------------------------------------------- (a) occlusion
ax = axes[0][0]
x = np.arange(len(order))
gly = [P[t]["glycan_clashes_3A"] for t in order]
nbr = [P[t]["neighbour_clashes_3A"] for t in order]
ax.bar(x - 0.2, np.maximum(gly, 0.4), width=0.38, color="#8c6d3f", label="glycan")
ax.bar(x + 0.2, np.maximum(nbr, 0.4), width=0.38, color="#5b6f8c",
       label="neighbouring protomer / gp41")
ax.set_yscale("log")
ax.set_ylim(0.35, 2000)
ax.set_xticks(x)
ax.set_xticklabels([SHORT[t] for t in order])
ax.set_ylabel("clashing atom pairs < 3 Å")
ax.legend(frameon=False, fontsize=7.5, loc="upper right", labelcolor=INK)
for i, t in enumerate(order):
    if gly[i] == 0 and nbr[i] == 0:
        ax.text(i, 0.45, "clear", ha="center", fontsize=6.5, color=MUTED)
dress(ax, "a", "Occlusion on the closed glycosylated trimer")

# ------------------------------------------------------- (b) epitope retention
ax = axes[0][1]
keep = [P[t]["epitope_sasa_retained_pct"] for t in order]
ax.bar(x, keep, width=0.62, color=[COLOR[t] for t in order])
ax.axhline(P["CD4"]["epitope_sasa_retained_pct"], color=MUTED, lw=0.9, ls=(0, (2, 2)))
ax.text(len(order) - 0.4, P["CD4"]["epitope_sasa_retained_pct"] + 1.5, "CD4 reference",
        ha="right", fontsize=7, color=MUTED)
ax.set_xticks(x)
ax.set_xticklabels([SHORT[t] for t in order])
ax.set_ylabel("epitope surface retained (%)")
ax.set_ylim(60, 104)
dress(ax, "b", "Epitope surface surviving assembly")

# ------------------------------------------------- (c) does the clash matter
ax = axes[1][0]
tags = [t for t in order if P[t]["ddG_fraction_in_glycan_clash"] is not None]
vals = [100 * P[t]["ddG_fraction_in_glycan_clash"] for t in tags]
xx = np.arange(len(tags))
ax.bar(xx, vals, width=0.58, color=[COLOR[t] for t in tags])
for i, v in enumerate(vals):
    ax.text(i, v + 2.5, "%.0f%%" % v, ha="center", fontsize=8, color=INK)
ax.set_xticks(xx)
ax.set_xticklabels([SHORT[t] for t in tags])
ax.set_ylabel("hotspot energy in glycan contact (%)")
ax.set_ylim(0, 112)
dress(ax, "c", "Does the clash reach the binding energy")

# ------------------------------------------------------- (d) state dependence
ax = axes[1][1]
st = sd["footprints"]
tags2 = [t for t in ("CD4", "17b", "D5", "D2", "D2p", "D3")
         if "local_fit_rmsd_A" in st[t]]
vals2 = [st[t]["local_fit_rmsd_A"] for t in tags2]
xx = np.arange(len(tags2))
ax.bar(xx, vals2, width=0.58, color=[COLOR[t] for t in tags2])
ax.axhline(sd["whole_protomer_fit_A"], color=MUTED, lw=0.9, ls=(0, (2, 2)))
ax.text(len(tags2) - 0.4, sd["whole_protomer_fit_A"] + 0.12,
        "whole protomer (%.1f Å)" % sd["whole_protomer_fit_A"], ha="right",
        fontsize=7, color=MUTED)
ax.set_xticks(xx)
ax.set_xticklabels([SHORT[t] for t in tags2])
ax.set_ylabel("footprint RMSD, closed vs open (Å)")
ax.set_ylim(0, max(vals2) * 1.35)
dress(ax, "d", "Epitope shape through CD4-induced opening")

fig.tight_layout(h_pad=2.6, w_pad=3.0)
out = os.path.join(FIG, "fig12_native_context.png")
fig.savefig(out, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote", os.path.basename(out))
