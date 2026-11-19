"""
Assemble the final figure set.

PyMOL renders each scene on a white field with whatever margin the camera
happens to leave. Those margins differ from panel to panel, so a figure built
from the raw renders looks ragged and wastes column width. Here every render is
cropped to its ink, placed on a grid, and given a panel letter and the labels
that belong in the image rather than in the caption.

Run after figures4.pml.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image

FIG = r"C:\Users\anant\Desktop\hiv_miniprotein_design\paper_analysis\04_figures"
D1, D2, D5, D4, D3 = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"
CD4F, ABF = "#99C2E0", "#F5CC8C"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})


def crop(name, pad=12):
    """Trim the white field around a render, leaving a small even margin."""
    im = Image.open(os.path.join(FIG, name)).convert("RGB")
    a = np.asarray(im)
    ink = (a < 247).any(axis=2)
    if not ink.any():
        return im
    ys, xs = np.where(ink)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, a.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, a.shape[1])
    return im.crop((x0, y0, x1, y1))


def panel(ax, name, letter, title=None, legend=None):
    ax.imshow(np.asarray(crop(name)))
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(-0.01, 1.0, letter, transform=ax.transAxes, fontsize=17,
            fontweight="bold", va="top", ha="right")
    if title:
        ax.set_title(title, fontsize=11.5, pad=6)
    if legend:
        ax.legend(handles=[Patch(facecolor=c, label=l) for l, c in legend],
                  loc="lower center", bbox_to_anchor=(0.5, -0.11),
                  ncol=min(len(legend), 4), frameon=False,
                  fontsize=9.5, handlelength=1.1, handleheight=1.1,
                  columnspacing=1.2)


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", os.path.basename(p))


# ---- Figure 1: the two receptor sites, and where the designs sit -----------
fig, axes = plt.subplots(1, 2, figsize=(11, 5.6))
panel(axes[0], "panel_sites.png", "a", "Receptor footprints on gp120",
      [("CD4 footprint", CD4F), ("17b / coreceptor footprint", ABF),
       ("CD4bs hotspots", "#B22222")])
panel(axes[1], "panel_designs_on_sites.png", "b", "D1 on the CD4 site, D5 at the coreceptor site",
      [("D1 (CD4-mimetic control)", D1), ("D5 (de novo)", D5)])
save(fig, "fig1_epitope_classes.png")

# ---- Figure 2: D1 recovers the CD4 binding mode ---------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
panel(axes[0], "panel_d1_cd4.png", "a", "D1 recovers the CD4 footprint (97%)",
      [("D1", D1), ("CD4 D1 domain", "#8C8C8C"), ("CD4 footprint", CD4F)])
panel(axes[1], "panel_hotspots.png", "b", "D1 engages the CD4bs hotspots",
      [("D1", D1), ("CD4bs hotspot sidechains", "#B22222")])
save(fig, "fig2_cd4_competition.png")

# ---- Figure 3: D5 blocks the coreceptor site ------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
panel(axes[0], "panel_d5_footprint.png", "a", "D5 beside the 17b footprint (covers 33%)",
      [("D5", D5), ("17b / coreceptor footprint", ABF)])
panel(axes[1], "panel_d5_17b.png", "b", "D5 collides with the 17b Fab (414 atom pairs < 3 Å)",
      [("D5", D5), ("17b Fab", "#8C8C8C")])
save(fig, "fig3_d5_17b_competition.png")

# ---- Figure 4: simultaneous engagement ------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 5.6))
panel(ax, "panel_simultaneous.png", "", "D1 and D5 occupy gp120 together",
      [("D1", D1), ("D5", D5), ("CD4 footprint", CD4F),
       ("17b footprint", ABF)])
save(fig, "fig4_simultaneous_binding.png")

# ---- Figure 5: the whole panel on one surface -----------------------------
fig, ax = plt.subplots(figsize=(7.6, 5.8))
panel(ax, "panel_landscape.png", "", "Epitope landscape of the design panel",
      [("D1", D1), ("D2", D2), ("D3", D3), ("D2'", D4), ("D5", D5)])
save(fig, "fig5_epitope_landscape.png")

# ---- Figure 9: tandem constructs ------------------------------------------
if all(os.path.exists(os.path.join(FIG, f))
       for f in ("panel_bispecific_D1_D5.png", "panel_bispecific_D1_D2.png")):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    panel(axes[0], "panel_bispecific_D1_D5.png", "a", "D1-linker-D5",
          [("Arm 1 (D1)", D1), ("Linker", "#8C8C8C"), ("Arm 2 (D5)", D5)])
    panel(axes[1], "panel_bispecific_D1_D2.png", "b", "D1-linker-D2",
          [("Arm 1 (D1)", D1), ("Linker", "#8C8C8C"), ("Arm 2 (D2)", D2)])
    save(fig, "fig9_bispecific_constructs.png")
