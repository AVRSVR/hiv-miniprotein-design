# Publication figure renders, v4.
#
# Changes from figures3: epitope footprints are painted onto the gp120 surface
# rather than shown only as a second protein, so each panel answers its question
# from the surface itself; framing is tight (panels are auto-cropped afterwards
# by compose_figures.py); and every scene uses one orientation of gp120 so the
# reader can carry the same mental frame from panel to panel.
#
# Okabe-Ito, CVD-safe: D1 #0072B2 | D2 #E69F00 | D5 #009E73 | D2' #D55E00 | D3 #CC79A7
# Run: pymol -cq figures4.pml

set ray_opaque_background, 1
bg_color white
set antialias, 2
set ray_shadows, 0
set cartoon_fancy_helices, 1
set surface_quality, 1
set specular, 0.12
set ambient, 0.34
set direct, 0.55
set depth_cue, 0
set ray_trace_fog, 0
set transparency_mode, 2

cd C:/Users/anant/Desktop/hiv_miniprotein_design/paper_analysis/04_figures

load scene/gp120_ref.pdb, gp120
load ../01_structures/1GC1.pdb, ref1gc1
load scene/binder_S310.pdb, D1
load scene/binder_Chai1_Best.pdb, D2
load scene/binder_Cyclic_Ring.pdb, D3
load scene/binder_model_4_chai.pdb, D4
load scene/binder_Latest_Best.pdb, D5

create CD4, ref1gc1 and chain C and polymer
create Ab17b, ref1gc1 and chain L+H and polymer
delete ref1gc1

set_color c_d1, [0.000, 0.447, 0.698]
set_color c_d2, [0.902, 0.624, 0.000]
set_color c_d5, [0.000, 0.620, 0.451]
set_color c_d4, [0.835, 0.369, 0.000]
set_color c_d3, [0.800, 0.475, 0.655]
set_color c_cd4foot, [0.60, 0.76, 0.88]
set_color c_17bfoot, [0.96, 0.80, 0.55]


run C:/Users/anant/Desktop/hiv_miniprotein_design/paper_analysis/scripts/face_helper.py

hide everything
color grey85, gp120
color grey55, CD4
color grey45, Ab17b
color c_d1, D1
color c_d2, D2
color c_d3, D3
color c_d4, D4
color c_d5, D5

# Receptor footprints, measured from the superposed reference complex rather
# than taken from the literature, so they match the competition analysis.
select cd4_foot, byres (gp120 within 5 of CD4)
select ab_foot,  byres (gp120 within 5 of Ab17b)
select hotspots, gp120 and resi 368+370+371+425+427+429+457+473+280
select CD4_D1, CD4 and resi 1-100

# One canonical view of gp120, reused by every surface panel.
hide everything
show surface, gp120
face cd4_foot or ab_foot, gp120, 3
view canon, store

# --------------------------------------------- Fig 1: the two receptor sites
hide everything
show surface, gp120
color grey85, gp120
color c_cd4foot, cd4_foot
color c_17bfoot, ab_foot
color firebrick, hotspots
view canon, recall
zoom gp120, 1
png panel_sites.png, width=2000, height=2000, dpi=300, ray=1

# --------------------------------------------- Fig 1b: where the designs sit
show cartoon, D1 D5
set cartoon_transparency, 0.0
view canon, recall
zoom gp120 or D1 or D5, 2
png panel_designs_on_sites.png, width=2000, height=2000, dpi=300, ray=1

# --------------------------------------------- Fig 2: D1 recovers the CD4 mode
hide everything
show surface, gp120
color grey85, gp120
color c_cd4foot, cd4_foot
show cartoon, CD4_D1
set cartoon_transparency, 0.55, CD4
show cartoon, D1
face cd4_foot, cd4_foot or D1, 4
png panel_d1_cd4.png, width=2000, height=1700, dpi=300, ray=1

# --------------------------------------------- Fig 3: D5 covers the 17b footprint
hide everything
show surface, gp120
color grey85, gp120
color c_17bfoot, ab_foot
show cartoon, D5
face ab_foot, ab_foot or D5, 2
png panel_d5_footprint.png, width=2000, height=1700, dpi=300, ray=1

# --------------------------------------------- Fig 3b: D5 against the 17b Fab
hide everything
show surface, gp120
set transparency, 0.45, gp120
color grey85, gp120
show cartoon, Ab17b
set cartoon_transparency, 0.5, Ab17b
show cartoon, D5
face ab_foot, D5 or ab_foot or Ab17b, 4
png panel_d5_17b.png, width=2000, height=1700, dpi=300, ray=1

# --------------------------------------------- Fig 4: both sites at once
hide everything
show surface, gp120
set transparency, 0.0, gp120
color grey85, gp120
color c_cd4foot, cd4_foot
color c_17bfoot, ab_foot
show cartoon, D1 D5
face cd4_foot or ab_foot, gp120 or D1 or D5, 3
png panel_simultaneous.png, width=2200, height=1800, dpi=300, ray=1

# --------------------------------------------- Fig 7: D1 on the CD4bs hotspots
hide everything
show cartoon, gp120
color grey80, gp120
show sticks, hotspots
color firebrick, hotspots
show cartoon, D1
orient hotspots
zoom hotspots, 8
png panel_hotspots.png, width=2000, height=1700, dpi=300, ray=1

# --------------------------------------------- Fig 5: whole panel on one surface
hide everything
show surface, gp120
set transparency, 0.15, gp120
color grey85, gp120
show cartoon, D1 D2 D3 D4 D5
view canon, recall
zoom gp120 or D1 or D2 or D3 or D4 or D5, 2
png panel_landscape.png, width=2200, height=1800, dpi=300, ray=1

print("FIGURES4 DONE")
