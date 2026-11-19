# Publication structural figures: designed binders on gp120 vs the CD4 footprint.
# Run: pymol -cq figures.pml

set ray_opaque_background, 1
bg_color white
set antialias, 2
set ray_shadows, 0
set cartoon_fancy_helices, 1
set cartoon_transparency, 0.0
set surface_quality, 1
set specular, 0.2
set ambient, 0.25
set direct, 0.6
set ray_trace_mode, 0
set depth_cue, 0

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

# CD4-binding-site hotspots on gp120 (HXB2)
select hotspots, gp120 and resi 368+370+371+425+427+429+457+473+280
color firebrick, hotspots

# ---------------------------------------------------------------- Figure 1
# Two distinct epitopes: CD4bs-blocking vs non-competitive
show cartoon, gp120
show surface, gp120
set transparency, 0.55, gp120
show cartoon, S310 BiSpec
set cartoon_transparency, 0.0, S310 BiSpec
show sticks, hotspots and not (name C+N+O)
orient gp120
turn y, -25
zoom gp120, 4
png fig1_two_epitopes.png, width=2400, height=1800, dpi=300, ray=1

# ---------------------------------------------------------------- Figure 2
# CD4 competition: S310 occupies the CD4 footprint, bispecific does not
hide everything
show cartoon, gp120
show surface, gp120
set transparency, 0.6, gp120
show cartoon, CD4
show cartoon, S310 BiSpec
orient gp120 or CD4
turn y, -20
zoom gp120 or CD4, 3
png fig2_cd4_competition.png, width=2400, height=1800, dpi=300, ray=1

# ---------------------------------------------------------------- Figure 3
# Simultaneous engagement: both designs bound at once, CD4 displaced
hide everything
show surface, gp120
set transparency, 0.35, gp120
show cartoon, S310 BiSpec
orient gp120
turn y, -25
turn x, 10
zoom gp120, 4
png fig3_simultaneous_binding.png, width=2400, height=1800, dpi=300, ray=1

# ---------------------------------------------------------------- Figure 4
# CD4bs hotspot engagement by S310 (close view)
hide everything
show cartoon, gp120
color grey70, gp120
show sticks, hotspots
color firebrick, hotspots
show cartoon, S310
show sticks, S310 and byres (S310 within 4.5 of hotspots)
orient hotspots or (S310 within 8 of hotspots)
zoom hotspots, 6
png fig4_hotspot_engagement.png, width=2400, height=1800, dpi=300, ray=1

# ---------------------------------------------------------------- Figure 5
# All designs mapped onto gp120: epitope landscape
hide everything
show surface, gp120
set transparency, 0.3, gp120
show cartoon, S310 BiSpec CycRing M4chai
orient gp120
turn y, -25
zoom gp120, 4
png fig5_epitope_landscape.png, width=2400, height=1800, dpi=300, ray=1

print "FIGURES DONE"
