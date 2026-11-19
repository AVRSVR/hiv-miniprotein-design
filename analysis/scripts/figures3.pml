# Final publication figures.
# Colours are the CVD-validated Okabe-Ito order used throughout the paper:
#   D1 #0072B2 blue | D2 #E69F00 orange | D5 #009E73 green
#   D4 #D55E00 vermillion | D3 #CC79A7 purple
# Run: pymol -cq figures3.pml

set ray_opaque_background, 1
bg_color white
set antialias, 2
set ray_shadows, 0
set cartoon_fancy_helices, 1
set surface_quality, 1
set specular, 0.15
set ambient, 0.32
set direct, 0.55
set depth_cue, 0

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

hide everything
color grey80, gp120
color grey60, CD4
color grey45, Ab17b
color c_d1, D1
color c_d2, D2
color c_d3, D3
color c_d4, D4
color c_d5, D5

select hotspots, gp120 and resi 368+370+371+425+427+429+457+473+280
color firebrick, hotspots
select CD4_D1, CD4 and resi 1-100

# --------------------------------------------- Fig 5: epitope landscape (all five)
hide everything
show surface, gp120
set transparency, 0.25, gp120
show cartoon, D1 D2 D3 D4 D5
orient gp120
turn y, -20
zoom gp120, 2
png fig5_epitope_landscape.png, width=2400, height=1800, dpi=300, ray=1

# --------------------------------------------- Fig 1: CD4bs vs non-CD4bs
hide everything
show surface, gp120
set transparency, 0.3, gp120
show cartoon, D1 D2 D5
show spheres, hotspots and name CA
set sphere_scale, 0.55, hotspots
orient gp120
turn y, -20
zoom gp120, 2
png fig1_epitope_classes.png, width=2400, height=1800, dpi=300, ray=1

# --------------------------------------------- Fig 2: CD4 competition
hide everything
show surface, gp120
set transparency, 0.5, gp120
show cartoon, CD4_D1
set cartoon_transparency, 0.4, CD4
show cartoon, D1 D2
orient gp120 or CD4_D1
turn y, -20
zoom gp120 or CD4_D1, 1.6
png fig2_cd4_competition.png, width=2400, height=1800, dpi=300, ray=1

# --------------------------------------------- Fig 2b: D1/CD4 overlap close-up
hide everything
show cartoon, gp120
color grey85, gp120
show sticks, hotspots
show cartoon, CD4_D1
set cartoon_transparency, 0.5, CD4
show cartoon, D1
orient D1 or CD4_D1
zoom D1 or CD4_D1, 2
png fig2b_d1_cd4_overlap.png, width=2400, height=1800, dpi=300, ray=1

# --------------------------------------------- Fig 3: D5 blocks the 17b/coreceptor site
hide everything
show surface, gp120
set transparency, 0.55, gp120
show cartoon, Ab17b
set cartoon_transparency, 0.45, Ab17b
show cartoon, D5
orient D5 or Ab17b
zoom D5 or Ab17b, 2.5
png fig3_d5_17b_competition.png, width=2400, height=1800, dpi=300, ray=1

# --------------------------------------------- Fig 4: simultaneous engagement D1 + D2
hide everything
show surface, gp120
set transparency, 0.3, gp120
show cartoon, D1 D2
orient gp120 or D1 or D2
turn y, -20
turn x, 8
zoom gp120 or D1 or D2, 2
png fig4_simultaneous_binding.png, width=2400, height=1800, dpi=300, ray=1

# --------------------------------------------- Fig 7: D1 hotspot engagement
hide everything
show cartoon, gp120
color grey75, gp120
show sticks, hotspots
color firebrick, hotspots
show cartoon, D1
orient hotspots
zoom hotspots, 7
png fig7_hotspot_engagement.png, width=2400, height=1800, dpi=300, ray=1

print "FIGURES3 DONE"
