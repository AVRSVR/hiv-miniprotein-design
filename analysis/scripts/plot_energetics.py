"""
Figures 7 and 8: interface hotspots and the escape-risk landscape.

Fig 7  per-design alanine-scan hotspot maps - which binder residues carry binding
Fig 8  escape-risk landscape - binder dependence on each gp120 contact plotted
       against how conserved that position is across HIV-1 Env
"""
import os, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
TBL = os.path.join(PROJECT, "paper_analysis", "05_tables")
FIG = os.path.join(PROJECT, "paper_analysis", "04_figures")

ORDER = ["Latest_Best", "S310", "Chai1_Best", "Cyclic_Ring"]
NAME = {"Latest_Best": "D5  de novo, coreceptor site (66 aa)",
        "S310": "D1  CD4-mimetic control (97 aa)",
        "Chai1_Best": "D2  de novo bispecific (73 aa)",
        "Cyclic_Ring": "D3  CD4-mimetic cyclic (13 aa)"}
COLOR = {"Latest_Best": "#009E73", "S310": "#0072B2",
         "Chai1_Best": "#E69F00", "Cyclic_Ring": "#CC79A7"}
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#dcdcdc"

# three-letter -> one-letter; resname[0] is wrong (ARG would render as "A")
AA3TO1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
          "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
          "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
          "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def one(resname):
    return AA3TO1.get(resname.upper(), "X")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.linewidth": 0.8, "figure.dpi": 300})

# ----------------------------------------------------------------- Figure 7
scans = {}
for tag in ORDER:
    p = os.path.join(TBL, f"alascan_{tag}.json")
    if os.path.exists(p):
        scans[tag] = json.load(open(p))

if scans:
    n = len(scans)
    fig, axes = plt.subplots(n, 1, figsize=(7.2, 2.3 * n))
    if n == 1:
        axes = [axes]
    for i, (ax, tag) in enumerate(zip(axes, [t for t in ORDER if t in scans])):
        d = scans[tag]
        res = sorted(d["residues"], key=lambda z: z["resid"])
        labels = [f"{one(z['resname'])}{z['resid']}" for z in res]
        vals = [z["ddG_kcal_mol"] for z in res]
        errs = [z.get("sem_kcal_mol", 0) or 0 for z in res]
        cols = [COLOR[tag] if v >= 1.0 else "#c9c9c9" for v in vals]
        x = np.arange(len(vals))
        ax.bar(x, vals, yerr=errs, color=cols, width=0.72,
               error_kw=dict(ecolor=MUTED, elinewidth=0.7, capsize=0))
        ax.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (2, 2)))
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
        ax.set_ylabel("ΔΔG (kcal/mol)")
        nh = sum(1 for v in vals if v >= 1.0)
        ax.text(-0.075, 1.12, "abcd"[i], transform=ax.transAxes, fontsize=16,
                fontweight="bold", va="top", ha="right")
        ax.set_title(f"{NAME[tag]}   —   {nh} hotspots ≥ 1 kcal/mol",
                     loc="left", fontsize=9.5, color=INK, pad=5)
        ax.grid(True, axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.margins(x=0.01)
    axes[-1].set_xlabel("Binder residue")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig7_alanine_scan.png"),
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote fig7_alanine_scan.png")

# ----------------------------------------------------------------- Figure 8
esc = {}
for tag in ORDER:
    p = os.path.join(TBL, f"escape_risk_{tag}.json")
    if os.path.exists(p):
        esc[tag] = json.load(open(p))

if esc:
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for tag, d in esc.items():
        pts = [c for c in d["contacts"]
               if c["conservation"] is not None and c["ddG_kcal_mol"] > 0]
        if not pts:
            continue
        xs = [100 * c["conservation"] for c in pts]
        ys = [c["ddG_kcal_mol"] for c in pts]
        ax.scatter(xs, ys, s=38, color=COLOR[tag], alpha=0.85,
                   edgecolor="white", linewidth=0.8, label=NAME[tag], zorder=3)
        # label the clearest liabilities
        for c in sorted(pts, key=lambda z: -(z["ddG_kcal_mol"] *
                                             (1 - z["conservation"])))[:3]:
            ax.annotate(f"{one(c['resname'])}{c['hxb2_pos']}",
                        (100 * c["conservation"], c["ddG_kcal_mol"]),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=7.5, color=INK)

    ax.axvspan(0, 80, color="#d9534f", alpha=0.055, zorder=0)
    ax.axvspan(95, 101, color="#5cb85c", alpha=0.07, zorder=0)
    ax.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (2, 2)), zorder=1)
    # keep the band label clear of the legend block in the upper centre
    ax.text(24, ax.get_ylim()[1] * 0.97, "variable positions\n(escape routes)",
            ha="center", va="top", fontsize=8, color="#b5423c")
    ax.text(97.5, ax.get_ylim()[1] * 0.96, "invariant\n(escape-resistant)",
            ha="center", va="top", fontsize=8, color="#3d8b40")
    ax.set_xlabel("Conservation of the gp120 position across HIV-1 Env (%)")
    ax.set_ylabel("Binder dependence, ΔΔG on removal (kcal/mol)")
    ax.set_title("Escape-risk landscape: heavy dependence on a variable position "
                 "is a cheap escape route",
                 loc="left", fontsize=10, color=INK, pad=8)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.56, 0.93), labelcolor=INK)
    ax.set_xlim(15, 103)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig8_escape_risk.png"),
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote fig8_escape_risk.png")
