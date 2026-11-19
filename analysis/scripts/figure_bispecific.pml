# Figure 8: modelled tandem bispecific construct bound to gp120.
# In the written complex, chain A is gp120 and chain B is the construct
# (arm 1 = residues 1-97, linker = 98-119, arm 2 = 120-end).

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
set cartoon_transparency, 0.0

cd C:/Users/anant/Desktop/hiv_miniprotein_design/paper_analysis/04_figures

set_color c_d1, [0.000, 0.447, 0.698]
set_color c_d2, [0.902, 0.624, 0.000]
set_color c_d5, [0.000, 0.620, 0.451]
set_color c_link, [0.45, 0.45, 0.45]

# ---------------- D1 -> D5 (attachment + coreceptor blockade)
load ../06_bispecific/S310__Latest_Best_with_gp120.pdb, bs1
hide everything
select gp1, bs1 and chain A
select arm1a, bs1 and chain B and resi 1-97
select link1, bs1 and chain B and resi 98-119
select arm1b, bs1 and chain B and resi 120-185

show surface, gp1
color grey80, gp1
set transparency, 0.35, gp1
show cartoon, arm1a link1 arm1b
color c_d1, arm1a
color c_link, link1
color c_d5, arm1b
orient bs1
turn y, -15
zoom bs1, 3
png panel_bispecific_D1_D5.png, width=2400, height=1800, dpi=300, ray=1
delete bs1

# ---------------- D1 -> D2
load ../06_bispecific/S310__Chai1_Best_with_gp120.pdb, bs2
hide everything
select gp2, bs2 and chain A
select arm2a, bs2 and chain B and resi 1-97
select link2, bs2 and chain B and resi 98-127
select arm2b, bs2 and chain B and resi 128-200

show surface, gp2
color grey80, gp2
set transparency, 0.35, gp2
show cartoon, arm2a link2 arm2b
color c_d1, arm2a
color c_link, link2
color c_d2, arm2b
orient bs2
turn y, -15
zoom bs2, 3
png panel_bispecific_D1_D2.png, width=2400, height=1800, dpi=300, ray=1

print "BISPECIFIC FIGURES DONE"
