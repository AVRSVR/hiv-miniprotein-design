# HIV gp120 Miniprotein Binder Design

Computational design and structural characterization of miniprotein binders against the HIV-1 envelope glycoprotein gp120, combining generative structure design, sequence design, structural validation, and molecular dynamics.

Epitopes were determined **from the structures themselves** rather than assumed from design intent — and the result was not what the campaign targeted. The designs engage surfaces *outside* the CD4-binding site, including one that blocks the coreceptor-binding site, and designs at independent surfaces can engage gp120 simultaneously.

## Approach

**1. Target preparation**
The HIV-1 envelope trimer structure (PDB [6V8X](https://www.rcsb.org/structure/6V8X)) was downloaded and processed: the protomer chain was isolated, cleaned and repaired (missing atoms/residues, protonation), and CD4-binding-site hotspot residues were mapped onto the cleaned structure to define the intended design target.

**2. Design generation**
Binder backbones and sequences were generated with **RFdiffusion** + **ProteinMPNN**, and independently with **Rosetta**-based protocols, followed by AlphaFold2 structural validation. Design–target complexes were produced by co-folding.

**3. Iterative refinement with Proteus**
Top candidates were carried through an iterative refinement loop using **[Proteus](https://proteus-5kpn.onrender.com/)** — our own in-house tool that automates the generate–check–refine cycle: it generates variants, runs a short MD stability check on each, and advances only those that pass.

**4. Structural characterization**
Each design's epitope was determined directly from its complex using two independent criteria (residue-contact mapping against the 14 canonical CD4bs hotspots, and a geometric test against the CD4bs centroid). Competition against host **CD4** and the coreceptor-site antibody **17b** was quantified by superposition onto the experimental gp120–CD4–17b ternary complex (PDB 1GC1) — no docking required.

**5. Molecular dynamics**
Every design was simulated for 2 ns (GBn2 implicit solvent, 310.15 K), tracking interface RMSD, contact persistence and binder fold stability, with MM/GBSA binding energies computed from the trajectories.

## Key results

**Designs bind outside the CD4-binding site.** The de novo designs contact **0 of 14** CD4bs hotspots, sitting 18–21 Å away, and show **zero** steric overlap with CD4.

**One design blocks the coreceptor site.** A 66-residue de novo binder overlaps the 17b Fab epitope (414 atom-level clashes) while showing no CD4 overlap — a distinct inhibition mechanism, acting *after* CD4 engagement, at a surface far less crowded than the heavily contested CD4bs. It is also the most conformationally stable binder in the panel (fold RMSD 1.35 Å, interface contacts at 101% of baseline over 2 ns).

**Designs at independent surfaces can bind simultaneously.** 7 of 10 design pairs show zero steric conflict, and tandem constructs joining two arms with a flexible linker were built and energy-relaxed without steric violations — a route to multi-epitope agents that raise the barrier to viral escape.

**Interfaces are of realistic scale.** All designs bury 1,254–2,589 Å², within the range spanned by the experimental gp120/CD4, gp120/b12 and gp120/17b interfaces (1,139–1,912 Å²).

![MD stability analysis](results/md_stability_analysis.png)

**MD trajectory (2 ns, water/ions hidden for visualization):** [media/md_2ns_run.mp4](media/md_2ns_run.mp4)

## Important caveats

These are computational results and hypothesis-generating only.

- **The complexes are predicted, not experimental.** Epitope assignments inherit any error in the co-folding models and require experimental validation (mutagenesis, competition ELISA, or structure determination).
- **Competition is assessed sterically** — necessary but not sufficient for functional competition; it captures neither affinity nor kinetics.
- **MD is short and single-replicate.** Run-to-run variance at 2 ns is substantial. MM/GBSA values carry no entropy term and support coarse ranking only — **no affinities should be inferred.**
- **Two scaffolds in the wider panel are CD4-mimetic** (derived from the human CD4 D1 domain) and serve as positive controls, not candidates — CD4-derived agents carry known liabilities including enhancement of infection in CD4-negative cells.
- **A single Env context.** One gp120 structure, one reference complex; cross-clade breadth is untested.
- **The bispecific's immune-engaging arm is not validated.** Only its gp120-facing interaction has been characterized.

## Next steps

- **Experimental validation** — recombinant expression, SPR/BLI against soluble gp120, and competition assays against soluble CD4 and 17b to test the predicted competition pattern directly.
- **Replicate and longer MD** to establish whether the stability ordering is robust.
- **Neutralization assays** against a multi-clade pseudovirus panel.
- **Immunogenicity and biocompatibility screening** before any further development.

## What's in this repo

```
src/            Core scripts: target download/prep, chain isolation, cleaning,
                hotspot mapping, and the explicit-solvent MD refinement loop
results/        MD stability plot and the designed structure (PDB)
media/          2 ns MD trajectory video (water/ions hidden)
docs/           Methods notes
```

## What's not included

This repo shows the pipeline and the headline results, not every intermediate artifact:
- The full set of trial designs and rejected candidates
- Intermediate scoring data and logs
- Proteus source code (internal tooling — see [proteus-5kpn.onrender.com](https://proteus-5kpn.onrender.com/))

## Requirements

- Python 3.9+
- BioPython, PDBFixer/OpenMM (target prep, MD)
- AlphaFold2 (structural validation)

> **Note on docking:** AutoDock Vina is unsuitable for binders of this size — it is a small-molecule engine, and a 70–100 residue binder has far too many rotatable bonds. Bound poses here come from co-folding, and CD4/17b competition is determined by superposition onto an experimental ternary complex, which requires no docking.
