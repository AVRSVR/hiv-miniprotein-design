# Publication structural figures (refined framing + labels).
# Run: pymol -cq figures2.pml

set ray_opaque_background, 1
bg_color white
set antialias, 2
set ray_shadows, 0
set cartoon_fancy_helices, 1
set surface_quality, 1
set specular, 0.15
set ambient, 0.30
set direct, 0.55
set depth_cue, 0
set label_size, 34
set label_color, black
set label_font_id, 7
set label_outline_color, white

cd C:/Users/anant/Desktop/hiv_miniprotein_design/paper_analysis/04_figures

load scene/gp120_ref.pdb, gp120
load scene/cd4_ref.pdb, CD4
load scene/binder_S310.pdb, S310
load scene/binder_Chai1_Best.pdb, BiSpec
load scene/binder_Cyclic_Ring.pdb, CycRing
load scene/binder_model_4_chai.pdb, M4chai

hide everything
color grey80, gp120
color forest, CD4
color marine, S310
color orange, BiSpec
color purple, CycRing
color salmon, M4chai
select hotspots, gp120 and resi 368+370+371+425+427+429+457+473+280
color firebrick, hotspots

# CD4 D1 is the gp120-binding domain; D2 projects away and only wastes frame
select CD4_D1, CD4 and resi 1-100

# ------------------------------------------------- Fig 2: CD4 competition
hide everything
show surface, gp120
set transparency, 0.55, gp120
show cartoon, CD4_D1
set cartoon_transparency, 0.45, CD4
show cartoon, S310 BiSpec
pseudoatom lab_cd4, pos=[0,0,0], label="CD4"
orient gp120 or CD4_D1
turn y, -20
zoom gp120 or CD4_D1, 1.5
png fig2_cd4_competition.png, width=2400, height=1800, dpi=300, ray=1

# ------------------------------------------------- Fig 2b: S310/CD4 overlap close-up
hide everything
show cartoon, gp120
color grey85, gp120
show sticks, hotspots
show cartoon, CD4_D1
set cartoon_transparency, 0.5, CD4
show cartoon, S310
orient S310 or CD4_D1
zoom S310 or CD4_D1, 2
png fig2b_s310_cd4_overlap.png, width=2400, height=1800, dpi=300, ray=1

# ------------------------------------------------- Fig 1: two epitopes
hide everything
show surface, gp120
set transparency, 0.4, gp120
show cartoon, S310 BiSpec
show spheres, hotspots and name CA
set sphere_scale, 0.6, hotspots
orient gp120
turn y, -20
zoom gp120, 2
png fig1_two_epitopes.png, width=2400, height=1800, dpi=300, ray=1

# ------------------------------------------------- Fig 3: simultaneous binding
hide everything
show surface, gp120
set transparency, 0.35, gp120
show cartoon, S310 BiSpec
orient gp120 or S310 or BiSpec
turn y, -20
turn x, 8
zoom gp120 or S310 or BiSpec, 2
png fig3_simultaneous_binding.png, width=2400, height=1800, dpi=300, ray=1

# ------------------------------------------------- Fig 5: epitope landscape
hide everything
show surface, gp120
set transparency, 0.3, gp120
show cartoon, S310 BiSpec CycRing M4chai
orient gp120
turn y, -20
zoom gp120, 2
png fig5_epitope_landscape.png, width=2400, height=1800, dpi=300, ray=1

print "FIGURES2 DONE"
