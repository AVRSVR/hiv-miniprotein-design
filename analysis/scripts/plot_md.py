"""
Figure 6: MD interface stability across the design panel.

Three stacked panels on a shared time axis (never a second y-scale):
  A  interface RMSD - does the binder stay in its bound pose
  B  retained interface contacts (% of initial) - does the interface survive
  C  binder fold RMSD - does the binder itself stay folded

Colours are the Okabe-Ito CVD-safe set in a validated order; line style carries
the same identity as a secondary encoding, so the series are separable in
greyscale and for colour-vision-deficient readers.
"""
import os, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
MD = os.path.join(PROJECT, "paper_analysis", "03_md")
FIG = os.path.join(PROJECT, "paper_analysis", "04_figures")
os.makedirs(FIG, exist_ok=True)

LABEL = {
    "S310":         "D1  CD4-mimetic control (97 aa)",
    "Chai1_Best":   "D2  de novo bispecific (73 aa)",
    "Cyclic_Ring":  "D3  CD4-mimetic cyclic (13 aa)",
    "model_4_chai": "D4  de novo (73 aa)",
    "Latest_Best":  "D5  de novo, coreceptor site (66 aa)",
}
ORDER = ["S310", "Chai1_Best", "Latest_Best", "model_4_chai", "Cyclic_Ring"]
# validated order: worst adjacent CVD dE 11.0
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"]
STYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]

INK = "#1a1a1a"
MUTED = "#6b6b6b"
GRID = "#dcdcdc"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.linewidth": 0.8, "figure.dpi": 300,
})

runs = []
for tag in ORDER:
    p = os.path.join(MD, f"{tag}_md_metrics.json")
    if not os.path.exists(p):
        continue
    m = json.load(open(p))
    fr = m.get("frames", [])
    if not fr:
        continue
    runs.append((tag, m, fr))

if not runs:
    raise SystemExit("no completed MD runs yet")

fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.2), sharex=True)
axA, axB, axC = axes

for i, (tag, m, fr) in enumerate(runs):
    t = np.array([x["time_ps"] for x in fr]) / 1000.0
    ir = np.array([x["iRMSD_A"] for x in fr])
    ct = np.array([x["contacts"] for x in fr], float)
    br = np.array([x["binder_RMSD_A"] for x in fr])
    c, ls = PALETTE[i % len(PALETTE)], STYLES[i % len(STYLES)]
    lab = LABEL.get(tag, tag)

    axA.plot(t, ir, color=c, ls=ls, lw=1.6, label=lab)
    axB.plot(t, 100.0 * ct / max(ct[0], 1), color=c, ls=ls, lw=1.6, label=lab)
    axC.plot(t, br, color=c, ls=ls, lw=1.6, label=lab)

axA.axhline(3.0, color=MUTED, lw=0.9, ls=(0, (2, 2)))
axA.text(0.995, 3.0, " 3 Å stability reference", transform=axA.get_yaxis_transform(),
         ha="right", va="bottom", fontsize=7.5, color=MUTED)
axA.set_ylabel("Interface RMSD (Å)")
axA.set_title("A   Binder retention in the bound pose", loc="left",
              fontsize=10, color=INK, pad=6)

axB.axhline(100.0, color=GRID, lw=0.9)
axB.set_ylabel("Interface contacts retained (%)")
axB.set_title("B   Interface persistence", loc="left", fontsize=10, color=INK, pad=6)

axC.set_ylabel("Binder fold RMSD (Å)")
axC.set_xlabel("Time (ns)")
axC.set_title("C   Binder conformational stability", loc="left",
              fontsize=10, color=INK, pad=6)

for ax in axes:
    ax.grid(True, color=GRID, lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.margins(x=0.01)

axA.legend(frameon=False, fontsize=7.8, loc="upper left",
           handlelength=2.6, labelcolor=INK)

fig.tight_layout()
out = os.path.join(FIG, "fig6_md_stability.png")
fig.savefig(out, bbox_inches="tight", facecolor="white")
print(f"wrote {out}  ({len(runs)} runs: {', '.join(t for t, _, _ in runs)})")
