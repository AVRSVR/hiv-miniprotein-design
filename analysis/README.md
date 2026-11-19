# Analysis pipeline

Every result in this directory is regenerated from the scripts in `scripts/`; nothing
is hand-entered. Each script writes its own JSON into `05_tables/`, the table builders
turn those into CSV and Markdown, and the plot scripts read the same JSON, so a number
never exists in two places.

## Design panel

| ID | Source file | Length | Origin | Epitope |
|----|-------------|--------|--------|---------|
| D1 | `God_Particle_S310` | 97 aa | CD4-mimetic (52.6% id to human CD4 D1) | CD4-binding site |
| D2 | `God_Particle_Chai1_Best` / `Strategy_E_Bispecific_Engager` | 73 aa | de novo | non-CD4bs, third site |
| D3 | `Strategy_D_Cyclic_Ring` | 13 aa | CD4-mimetic (100% id over 11 aa) | CD4-binding site (partial) |
| D4 | `model_4_chai` | 73 aa | de novo | CD4bs-peripheral |
| D5 | `Latest_Best` | 66 aa | de novo | **CD4i / coreceptor site** |

D1 and D3 are CD4-derived and serve as positive controls: they establish that the
analysis detects CD4bs engagement when it is genuinely present. D2, D4 and D5 are
de novo and carry the claims under test.

`docked_complex` is analysed for interface metrics and MD only. Its gp120 frame
superposes onto the 1GC1 reference at 26 Å RMSD, so its epitope cannot be assigned
and it is excluded from the competition and co-binding results.

## What is here

`scripts/` is one script per result, `05_tables/` the derived tables as CSV and
Markdown beside the JSON they were built from, `04_figures/` every figure, and
`01_structures/` the complexes and reference structures each measurement reads.

Headline results:
- **D5 (de novo, 66 aa) occludes the CD4i/coreceptor site** — 414 atom pairs within
  3 A of the 17b Fab, zero with CD4 — at a surface far less contested than the CD4bs.
  The occlusion is steric rather than epitope-centred: D5 covers 33% of the 17b
  footprint (5/15 residues) and sits mostly on adjacent surface. The CD4-mimetic
  control D1, by contrast, covers 97% of the CD4 footprint, which is what validates
  the measurement.
- D5 is the best candidate on every axis assessed: most conformationally stable
  (fold RMSD 1.35 A, contacts 101% of baseline), 18 alanine-scan hotspots totalling
  72.6 kcal/mol distributed across the interface, and the best developability
  profile (lowest exposed hydrophobic fraction 0.156, no glycosylation sequons,
  no unpaired cysteines).
- **Escape-risk map** names specific liabilities: D5 leans on Arg444 (+6.6 kcal/mol
  at a position only 40% conserved). Escape-resistant contact fraction ranks
  D1 (26%) > D5 (18%) > D2 (0%).
- 7 of 10 design pairs occupy gp120 simultaneously; two tandem constructs built and
  relaxed with zero steric violations and arms displaced only 0.32-0.42 A.
- **Redesign round is a negative result.** Four point mutants aimed at D5's Arg444
  liability all *increased* dependence on the variable positions (12.9 kcal/mol for
  the parent vs 14.0-31.6 for the mutants). The liability is positional, not
  sidechain-level, so relieving it needs backbone redesign.
- **The interface is van der Waals-driven, not a charge smear.** D5's electrostatic
  share of net binding (-0.76) sits inside the range set by the experimental
  CD4-gp120 (-0.87) and b12-gp120 (-0.53) interfaces scored identically. This
  answers the obvious objection to an E/K-rich helical hairpin; it does not
  establish specificity, which needs a counter-screen.
- **Run-to-run spread is ~4x the within-run SEM.** Three D5 trajectories give
  -75.6 +/- 3.1 kcal/mol against a per-trajectory SEM of 0.7. Differences below
  ~3 kcal/mol between designs are not real.

Two corrections worth carrying forward:
- D1 and D3 are **CD4-mimetic** (52.6% and 100% identity to human CD4) — positive
  controls, not candidates.
- **D2 and D2' are the same 73-residue sequence** in two predicted poses 11.3 A
  apart, not two designs. Retained deliberately as a pose-uncertainty control.
- D3 is not supported dynamically: it leaves its pose within 100 ps.

## Pipeline

| Step | Script | Output |
|------|--------|--------|
| Inventory all candidate structures | `inventory.py` | `00_inventory/inventory.json` |
| Split complexes, characterize interfaces | `interface_analysis.py` | `01_structures/interface_metrics.json` |
| Independent epitope check (geometry) | `verify_epitope.py` | `01_structures/epitope_verification.json` |
| CD4 / 17b competition + co-binding | `cd4_competition.py` | `05_tables/cd4_competition.json` |
| Place all binders in the 1GC1 frame | `export_superposed.py` | `04_figures/scene/*.pdb` |
| Interface-stability MD (2 ns each) | `run_md.py`, `md_queue.py` | `03_md/*_md_metrics.json` |
| MM/GBSA binding energies | `mmgbsa.py` | `05_tables/mmgbsa.json` |
| Interface benchmark vs experimental | `benchmark_interfaces.py` | `05_tables/interface_benchmark.json` |
| Cross-strain epitope conservation | `conservation.py` | `05_tables/epitope_conservation.json` |
| Tandem construct geometry | `bispecific_design.py` | `05_tables/bispecific_constructs.json` |
| Build + relax 3D bispecific models | `build_bispecific.py` | `06_bispecific/*_construct.pdb` |
| Computational alanine scan (binder) | `alanine_scan.py` | `05_tables/alascan_*.json` |
| Escape-risk map (gp120 side + conservation) | `escape_risk.py` | `05_tables/escape_risk_*.json` |
| Developability profiling | `developability.py` | `05_tables/developability.json` |
| Replicate MD | `replicate_queue.py` | `03_md/*_rep*_md_metrics.json` |
| D5 redesign round | `redesign_d5.py` | `07_redesign/redesign_d5.json` |
| Rescore redesign from saved trajectories | `rescore_d5.py` | `07_redesign/redesign_d5.json` |
| Tables 10-11 | `build_table10.py` | `05_tables/table1[01]_*.csv` |
| Tables 7-9 | `build_tables2.py` | `05_tables/table[789]_*.csv` |
| Figures 9-10 | `plot_energetics.py` | `04_figures/fig9,fig10*.png` |
| Publication tables | `build_tables.py` | `05_tables/table*.csv`, `.md` |
| Structural renders | `figures4.pml`, `face_helper.py` | `04_figures/panel_*.png` |
| Compose labelled multi-panel figures | `compose_figures.py` | `04_figures/fig[1-5,9]*.png` |
| Footprint overlap vs CD4/17b epitopes | `footprint_overlap.py` | `05_tables/footprint_overlap.json` |
| Redesign figure | `plot_redesign.py` | `04_figures/fig10_redesign.png` |
| Interface energy decomposition | `energy_decomposition.py` | `05_tables/energy_decomposition.json` |
| Glycan sequon proximity per footprint | `glycan_occlusion.py` | `05_tables/glycan_occlusion.json` |
| Per-position conservation (cached) | `conservation_core.py` | `05_tables/conservation_by_position.json` |
| Re-attach conservation to escape scans | `rescore_escape_risk.py` | `05_tables/escape_risk_*.json` |
| Tandem MD, per-arm | `analyze_tandem.py` | `05_tables/tandem_md.json` |
| Decomposition figure + Table 13 | `plot_decomposition.py` | `04_figures/fig11*.png` |
| MD figure | `plot_md.py` | `04_figures/fig6_md_stability.png` |

Run order: `inventory` → `interface_analysis` → `verify_epitope` → `cd4_competition`
→ `export_superposed` → figures/tables. MD is independent and can run any time after
`interface_analysis` (which produces the split PDBs).

## Method notes

**Why no docking.** AutoDock Vina is a small-molecule engine; the binders here have
far too many rotatable bonds for it, which is why the original `phase9` docking wrote
0.0 for every row. Bound poses come from co-folding instead, and CD4 competition is
determined by superposition onto the experimental gp120–CD4–17b complex (PDB 1GC1),
which needs no docking at all.

**Why superposition is fitted locally.** Our complexes are full-length gp120 while
1GC1 is the core construct. Fitting globally drags loop differences into the region
the competition question depends on, so the fit is restricted to the 37 gp120 residues
within 5 Å of CD4. Complexes fitting worse than 10 Å RMSD are rejected.

**Why MD restrains distal gp120.** gp120 is natively trimer-embedded and glycosylated,
so its distal scaffold is not free to move in reality. Restraining it prevents
implicit-solvent unfolding artifacts while leaving the interface unrestrained.

**Trajectory output must disable periodic-box enforcement.** These systems are
non-periodic (implicit solvent, NoCutoff), but OpenMM still carries a default 2 nm
box and `DCDReporter` wraps coordinates into it by default. The complexes are ~10 nm
across, so the saved trajectory is torn apart while the simulation itself is
unaffected — in-memory metrics stay correct, but anything computed from the DCD
(MM/GBSA, trajectory figures) is garbage. `enforcePeriodicBox=False` is required.

**DCD handles must be closed before the trajectory is read back.** On Windows the
`DCDReporter` keeps its file open after the run, and mdtraj then fails with
"premature end of file". The first redesign pass lost every result this way after
the MD itself had completed correctly; `rescore_d5.py` exists to recover such runs
without repeating the simulation.

**PDBFixer renumbers chains.** `applyMutations` and the repair steps renumber each
chain sequentially from 1, so HXB2 gp120 numbering no longer addresses the intended
residue in the rebuilt topology. Any scan carried over from `escape_risk.py` must
remap positionally first; without it the scan silently reports the wrong residues.

**Clash counts and footprint overlap answer different questions.** "Would the
receptor be blocked?" is a clash count; "does the binder sit on the receptor's
epitope?" is a footprint overlap. D5 scores high on the first and low on the second,
so only the blocking claim is made. `footprint_overlap.py` exists to keep the two
apart.

**`orient` is the wrong camera for a surface patch.** It aligns principal axes,
which leaves a patch on a roughly spherical protein edge-on with anything bound to
it hidden behind the bulk. `face_helper.py` places the camera on the ray from the
molecule's centre of mass through the patch's, which is what the epitope figures
need.

**MM/GBSA here ranks, it does not measure.** Single-trajectory, no entropy term,
implicit solvent. Useful for coarse ordering of related complexes; not affinities.

**Epitope assignment is measured, not assumed.** Two independent criteria are used —
residue contact mapping against the 14 CD4bs hotspots, and a geometric test against
the CD4bs centroid computed from each complex's own hotspot coordinates. They agree
on every design. Fold-confidence metrics report whether a sequence folds as intended,
not whether it engages the intended surface, so they are not used for this purpose.

## Known data issues in the upstream project

Recorded here so they are not silently inherited:

- The `blast_screening_results.csv` human off-target screen passes all 520 designs,
  including D1 (52.6% identical to human CD4) and D3 (100% over 11 residues). The
  screen does not detect real human-protein similarity and should be re-run.
- `iedb_results.csv` contains API `ERROR` rows; MHC-I immunogenicity coverage is
  therefore incomplete.
- `af2_rosetta_scores.csv` and `af2_rfdiffusion_*_scores.csv` are not comparable:
  the Rosetta arm was scored as single-chain monomers (pLDDT 0–100, intra-chain PAE)
  while the RFdiffusion arm was scored as two-chain complexes (pLDDT 0–1, `i_pae`
  inter-chain). No cross-arm comparison is drawn from it.
- `phase9_docking_results.csv` contains 0.0 for every binding energy (see "Why no
  docking"). Not used.

### Glycans are stripped, so sequon proximity has to be measured

Preparation removes heteroatoms, which on gp120 deletes every N-glycan. That makes
the MD tractable and every interface here is computed on bare surface, so a
footprint sitting where a glycan would sit looks perfectly good on every protein
metric. `glycan_occlusion.py` finds N-X-S/T sequons from the model's own sequence
(requiring consecutive residue numbers, so a sequon is never called across a chain
break) and profiles each footprint against the CD4 and 17b footprints, which bind
glycosylated gp120 and therefore define the tolerable range. D2 draws 58% of its
load-bearing interface energy from the N295 and N332 sequons, so most of its
interface is an artifact of the deletion; D5 draws 0%. Neither natural partner
contacts a sequon asparagine at all.

### A hand-built linker will not survive 4 fs

The tandem construct's poly-glycine linker is placed geometrically, and four of its
CA-CA distances come out near 3.1 A rather than 3.8 A - compressed, cis-like peptide
geometry. Bond lengths and nonbonded contacts are fine, so the problem is invisible
to a clash check, but it integrates to NaN. `run_md2.py` takes `MD_MIN_ITERS` and
`MD_RAMP_PS` from the environment for a longer minimisation and staged heating;
defaults reproduce the published panel runs unchanged. `MD_DT_FS` and `MD_HMASS`
are overridable for the same reason; the tandem ran at 2 fs with unrepartitioned
masses. Minimisation converges cleanly (-4.0e4 kJ/mol, largest residual forces about
1030 kJ/mol/nm on ordinary hydroxyl rotamers), so this is integration stability, not a
broken topology. The cis-like distances are still there after 2 ns, as expected -
isomerisation is far slower than the trajectory - so the linker conformation in these
models is not interpretable and needs rebuilding before anything is expressed.

Truncating around a 185-residue binder also strands isolated receptor residues, which
capping turns into free zwitterionic amino acids with no amber template. `run_md2.py`
now drops fragments under three residues and logs each one with its distance to the
binder; for the tandem that was 7 residues, all 6.6-10.5 A away and so outside contact
range.

### The tandem's arms hold, its linker does not behave

Restrained minimisation cannot show a two-armed construct works, because the restraints
hold the arms in place. Under 2 ns of free dynamics both arms stay engaged (D1 93%
contact retention at 2.36 A interface RMSD, D5 104% at 3.22 A). `analyze_tandem.py`
scores each arm separately against its own nearest receptor residues, because a
whole-construct iRMSD can look excellent while one arm releases and the other rotates
to compensate. Contact retention is quoted within the run: absolute counts depend on
how gp120 was truncated and do not compare across runs.

The linker is the part that misbehaves. The built model has it 10.4 A clear of gp120;
within equilibration it adsorbs onto the surface and stays at 2.87 A with about 70
heavy-atom contacts. Not a clash, and it detached neither arm, but the static model
gave no sign of it.

### Conservation came from 47 sequences and that was not enough

The escape argument originally rested on 47 reviewed UniProt Env entries, because
reviewed coverage of HIV-1 Env is tiny. Expanding to every full-length group M record
(4,497 fetched, 4,284 retained) moved the escape-resistant fractions enough to swap
two designs' places: D5 fell from 18% to 12% and D2 rose from 0% to 20%, so the lead
design is now the least escape-resistant of the three. Per-contact liabilities barely
moved (Arg444 40% -> 36%), so the specific redesign targets were sound while the
summary fractions were not.

Deposition is uneven - the largest set of near-identical records held 103 sequences -
so pooling over-weights whatever was sequenced most. Subtype balancing is impossible
here: 4,439 of 4,497 records carry no subtype in their metadata. `conservation_core.py`
instead greedy-clusters at 98% identity over the aligned HXB2 columns and scores one
representative per cluster. It also now owns the alignment outright; `conservation.py`
and `escape_risk.py` each used to carry their own copy, which was two chances to drift.
A LANL filtered web alignment is still the better source.

Only conservation depends on the sequence set, so `rescore_escape_risk.py` re-attaches
it to the stored scans rather than re-running the decoupling calculation, whose
per-contact dG values cannot change.

### The target came from a trimer, so the trimer went back

Every interface here was computed on one deglycosylated gp120 protomer. The protomer
came from 6V8X, and four of the five co-folded complexes still carry it at its
deposited coordinates (chain A superposes on 6V8X chain A at 0.00 A), so
`native_context.py` restores the neighbouring protomers, gp41, the 105 modelled sugar
residues and the bound VRC01 Fab with no superposition at all. The references validate
the test in both directions: CD4, which binds the closed trimer, makes one clash with
the quaternary structure, while 17b, which cannot bind it until CD4 opens it, makes 774.

D5 scores with CD4 (three clashes) and its conflict is with glycans at N262 and N295,
touching 23% of its hotspot energy. D2 draws 96% of its hotspot energy from
glycan-contacting residues, which is the same verdict `glycan_occlusion.py` reached from
sequon distances on a stripped model, by a method that never looked at a sugar atom.

D1 is the exception and finding that out mattered: its complex carries a CD4-bound
gp120 (16.5 A from 6V8X chain A over 434 CA) while the other four use the closed
protomer, so the panel is not conformationally uniform. D1, CD4 and 17b are placed by
footprint fit and their counts inherit that error.

### Two states of one isolate, so the comparison means something

`state_dependence.py` compares 6V8X chain A against 6CM3 chain D, both BG505, so no
isolate or construct difference enters. CD4's and 17b's footprints rearrange by 6.64 and
4.90 A through CD4-induced opening, which is what it means for an epitope to be formed by
the transition; D5's rearranges by 2.52 A. Carrying each binder onto the open protomer
also re-tests the central competition claim in the native assembly rather than by
superposing the 1GC1 core: D5 clashes with the 17b Fv and not with CD4.

### Rigid-body search bounds the glycan problem without settling it

`glycan_accommodation.py` searches 40,000 small rigid displacements for a pose that
clears every sugar. None exists within 4 A for any design that starts in conflict. Both
binder and glycans are rigid in that search and neither is rigid in reality, so it bounds
what nudging can fix rather than predicting that the binder fails. The experiment that
settles it is binding against glycosylated and deglycosylated gp120.

### The designs do not land on the one druggable pocket already known

`temsavir_site.py` carries BMS-626529 from 5U7O into the design frame by superposing its
gp120 onto 6V8X chain A (1.27 A over 429 CA, same isolate). The pocket is the canonical
24 residues. D5 is 13.6 A away and shares none of them, D2 18.3 A, D2' 13.1 A. D3, the
cyclized CD4 loop, closes to 6.8 A and shares three, which is what a CD4-mimetic should
do since the Phe43 cavity overlaps the drug site. So the de novo designs are not an
expensive route to a site medicinal chemistry already has, and nothing structural stops
a D5-class agent being combined with fostemsavir.

### A counter-screen, at last

`counter_screen.py` places D5 face-on at random surface points of human serum albumin
(1AO6), hen egg lysozyme (1DPX), ubiquitin (1UBQ) and gp120 itself, minimises every
placement for the same 1,500 steps and scores it with the same MM/GBSA. The designed
pose scores -50.9 kcal/mol; the best of eight placements is -11.1 on gp120, -26.2 on
albumin, -19.3 on lysozyme and -14.7 on ubiquitin, with medians of -5 to -11 everywhere.

So the designed interface is in a different class from arbitrary contact. The honest
other half: arbitrary placements on albumin score better than arbitrary placements on
gp120, which is unsurprising for a large groove-rich carrier protein against a flat
viral surface, but it means the screen bounds promiscuity rather than demonstrating
preference. Eight placements is not a docking search; every off-target number is an
upper bound.

The face-on constraint matters. Approaching from a random direction in a random
orientation makes first contact on one protruding atom almost every time, which yields
glancing touches rather than placements a binder could adopt. Presenting the binder's
own contacting face at the surface normal is what makes the decoys comparable.

### Stop using the surrogate

`ccr5_competition.py` carries CCR5 out of 6MEO, the gp120-CD4-CCR5 complex, into the
CD4-bound open trimer, so the design can be measured against the receptor instead of
against 17b.

Two traps had to be cleared first. 6MEO's gp120 is not HXB2-numbered and the offset
drifts along the chain, so pairing residues by number compares different residues to
each other: only 28% of shared numbers carry the same amino acid. Matching by sequence
alignment instead gives 365 pairs at 81% identity. And the two gp120s are not one rigid
body. Over all pairs they fit at 4.60 Å, with the inner domain 3.60 Å out and V1V2
4.18 Å, while the bridging sheet, outer domain, V3 and beta20-23 all agree to about an
angstrom. The parts that disagree are the parts that move on opening; the parts that
agree are the ones the receptor reads. So the transfer is local, over the 66 gp120
residues within 10 A of CCR5, fitting at 1.12 A.

The carry is checked against the interface it has to reproduce rather than asserted.
Natively CCR5 makes 504 contacts and 15 clashes with its own gp120; carried, it makes
235 and 14. Clashes preserved, contacts halved, so the pose transferred and the
interface loosened by about an angstrom. Read the contact numbers as relative.

What lands is the coreceptor site as the literature describes it, which nothing in the
procedure forced: the beta2-beta3 stem, the V3 base and crown, the beta20-beta21
strands, with sulfotyrosines Tys10 and Tys14 touching gp120.

Controls in opposite directions. CD4, which binds gp120 alongside CCR5 in the same
structure, makes no clash. 17b, which blocks coreceptor engagement, makes 357, with 248
against the N-terminal segment. D5 makes 47, with 20 against the N-terminal segment, and
covers 28.6% of the receptor's footprint against 17b's 32.1%.

The unflattering half: only 3 of D5's 17 load-bearing contacts are residues CCR5 itself
reads. It blocks by bulk, not by competing for what the receptor needs, so escape at
most of its contacts would not cost the virus its coreceptor.

### Conservation is not coverage

`strain_coverage.py` asks a question per-position conservation cannot answer: how much
of each interface a given circulating virus actually presents. `conservation_core.py`
now caches the clustered alignment as a matrix as well as a summary, so the 2,948
representatives can be scored one at a time.

Two distinctions drive it. Conservation counts the modal residue; a design holds the
target's residue, and those differ more often than one expects. And a binder needs its
contacts intact in the same virus, which multiplying per-position frequencies would
get wrong under linked selection.

Weighting each contact by alanine-scan ddG, a median virus presents 57.0% of D1's
interface, 61.0% of D5's and 65.6% of D2's, against 82.9% for CD4's epitope, 88.2% for
17b's and 89.3% for CCR5's footprint. The null is 2,000 contiguous patches of exposed
surface of matched size, which give 76.5% (5th to 95th percentile 41.2 to 94.1). D5 sits
at the 23rd percentile of that; the CCR5 footprint at the 93rd.

The shortfall is localised, not diffuse: five of D5's seventeen contacts carry 62% of it.
Two are the liabilities already known from conservation, Arg444 and Lys305. Two more are
a different kind of problem the conservation score cannot see, because at Gln440 and
Val442 the design target carries a residue found in 7% of the epidemic while the modal
residue is something else. Nothing has to mutate for those to cost binding.

An intuition that turned out false, recorded because it is the obvious one: the designs
do rest more heavily on the variable loops than the natural partners (42-63% of weight
against 12-17%), but that is not what costs them coverage. D5's ten loop contacts retain
74.7% and its seven non-loop contacts retain 35.4%. The damage is in the 440-446 stretch,
not in V1V2 or V3.

One generic trap worth repeating: the two stabilized BG505 constructs here differ at 4
of 363 shared positions, and at 204 the design target carries a valine where the other
construct and 98% of the epidemic carry alanine. Scoring 17b's epitope against the wrong
construct moved its coverage by 6 points. Every contact is now scored against the
residue in the structure where that contact was measured.
