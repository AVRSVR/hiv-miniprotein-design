"""
Tables 21-23: coreceptor competition and cross-strain coverage.

Run: python build_tables4.py
"""
import io
import json
import os

import pandas as pd

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(PA, "05_tables")


def load(name):
    with io.open(os.path.join(T, name), encoding="utf-8") as fh:
        return json.load(fh)


def write(df, stem, note):
    df.to_csv(os.path.join(T, stem + ".csv"), index=False)
    with io.open(os.path.join(T, stem + ".md"), "w", encoding="utf-8") as fh:
        fh.write(df.to_markdown(index=False))
        fh.write("\n\n" + note + "\n")
    print("\n" + stem)
    print(df.to_string(index=False))


# ------------------------------------------------------------------ Table 21: CCR5
d = load("ccr5_competition.json")
LAB = {"CD4": "CD4 (must not clash)", "17b": "17b Fab (must clash)",
       "D5": "D5 (de novo, coreceptor site)", "D2": "D2 (de novo bispecific)",
       "D2p": "D2′ (alternative pose)", "D3": "D3 (CD4-mimetic cyclic)"}
rows = []
for tag in ("CD4", "17b", "D5", "D2", "D2p", "D3"):
    r = d["partners"][tag]
    rows.append({"Partner": LAB[tag],
                 "gp120 footprint": r["footprint_residues"],
                 "Clashes < 3 Å": r["ccr5_clashes_3A"],
                 "Clashes < 2.5 Å": r["ccr5_severe_2.5A"],
                 "Closest approach (Å)": r["ccr5_min_dist_A"],
                 "Contacts < 4.5 Å": r["contacts_with_ccr5_4.5A"],
                 "Clashes with N-terminal segment": r["nterm_clashes_3A"],
                 "Clashes with receptor body": r["body_clashes_3A"],
                 "CCR5 footprint covered (%)": r["pct_of_ccr5_footprint_covered"]})
v = d["carry_validation"]
write(pd.DataFrame(rows), "table21_ccr5_competition",
      "CCR5 taken from %s and carried onto %s by a local superposition over the %d gp120 "
      "residues within %.0f Å of the receptor, fitting at %.2f Å; the whole-domain fit over "
      "%d sequence-matched pairs is %.2f Å, with the error concentrated outside the "
      "receptor's footprint. The transfer preserves the interface it came from: %d "
      "contacts and %d clashes natively against %d and %d after the carry. The receptor's "
      "gp120 footprint is %d residues. Sulfotyrosines Tys10 and Tys14 contact gp120 in the "
      "carried pose. The designs are carried onto the open protomer by their own footprint "
      "fits as in Table 17, so every number inherits both transforms, and all molecules "
      "are rigid."
      % (d["coreceptor_source"], d["open_trimer"], d["local_fit_residues"],
         d["local_fit_shell_A"], d["local_fit_A"], d["matched_pairs"], d["global_fit_A"],
         v["native_contacts_4.5A"], v["native_clashes_3A"], v["carried_contacts_4.5A"],
         v["carried_clashes_3A"], len(d["ccr5_footprint_on_gp120"])))

# -------------------------------------------------------- Table 22: strain coverage
c = load("strain_coverage.json")
ORDER = [("D1", "D1 hotspots (ΔΔG weighted)"), ("D5", "D5 hotspots (ΔΔG weighted)"),
         ("D2", "D2 hotspots (ΔΔG weighted)"),
         ("CD4 epitope", "CD4 epitope (equal weights)"),
         ("17b epitope", "17b epitope (equal weights)"),
         ("CCR5", "CCR5 footprint (equal weights)")]
rows = []
for key, lab in ORDER:
    r = c["sets"].get(key)
    if not r:
        continue
    rows.append({"Interface": lab,
                 "Positions": r["positions_scored"],
                 "Median retained (%)": r["median_retained_pct"],
                 "Mean retained (%)": r["mean_retained_pct"],
                 "Strains ≥ 90% (%)": r["strains_retaining_90pct_or_more"],
                 "Strains ≥ 75% (%)": r["strains_retaining_75pct_or_more"],
                 "Strains ≥ 50% (%)": r["strains_retaining_50pct_or_more"],
                 "Positions where target is not modal":
                     len(r["positions_where_bg505_is_not_modal"]),
                 "Weight in V1V2/V3/V4/V5 (%)": r["weight_in_variable_loops_pct"]})
n5 = c["null_model"].get("D5", {})
n28 = c["null_model"].get("CCR5", {})
write(pd.DataFrame(rows), "table22_strain_coverage",
      "Coverage over %d cluster representatives drawn from %d sequences at %.0f%% "
      "identity. A position counts as retained when the sequence carries the same residue "
      "as the structure in which that contact was measured: 6V8X chain A for the designs, "
      "6CM3 chain D for the reference footprints, which differ at 4 positions. Design "
      "contacts are weighted by alanine-scan ΔΔG above %.1f kcal/mol; reference footprints "
      "are weighted equally. Alignment gaps count as mismatches. The null is %d "
      "contiguous surface patches of matched size: %d residues give a median of %.1f%% "
      "(5th to 95th percentile %.1f to %.1f), %d residues %.1f%% (%.1f to %.1f). D5's set "
      "sits at the %.1fth percentile of its null and the CCR5 footprint at the %.1fth. A "
      "matching residue is necessary, not sufficient, so these are upper bounds on breadth "
      "and not neutralization predictions. The last column is the share of each "
      "interface that rests on the variable loops, which is where the "
      "coverage difference comes from."
      % (c["representatives"], c["sequences_behind_them"],
         100 * c["cluster_identity"], c["hotspot_threshold_kcal"], n5.get("draws", 0),
         n5.get("size", 0), n5.get("median_of_medians_pct", 0), n5.get("p5_pct", 0),
         n5.get("p95_pct", 0), n28.get("size", 0), n28.get("median_of_medians_pct", 0),
         n28.get("p5_pct", 0), n28.get("p95_pct", 0),
         n5.get("percentile_of_observed", 0), n28.get("percentile_of_observed", 0)))

# ------------------------------------------------- Table 23: position detail for D5
r = c["sets"]["D5"]
rows = []
for q in sorted(r["per_position"], key=lambda x: -x["weight_kcal"]):
    rows.append({"Position": q["position"],
                 "Target residue": q["bg505_residue"],
                 "ΔΔG (kcal/mol)": q["weight_kcal"],
                 "Target residue frequency (%)": round(100 * q["bg505_frequency"], 1),
                 "Modal residue": q["modal_residue"],
                 "Modal frequency (%)": round(100 * q["modal_frequency"], 1),
                 "Target is modal": "yes" if q["bg505_is_modal"] else "no",
                 "In a variable loop": "yes" if q["in_variable_loop"] else "no",
                 "Share of shortfall (%)": q["share_of_shortfall_pct"]})
write(pd.DataFrame(rows), "table23_d5_position_coverage",
      "Every load-bearing contact of D5, ordered by ΔΔG. The two frequency columns are "
      "the quantity that matters for a design and the quantity conservation reports: how "
      "often the epidemic presents the residue the design is touching, and how often it "
      "presents the most common one. Where they differ, the position can be well conserved "
      "and the design still mismatched. Frequencies inside V1V2, V3, V4 and V5 rest on a "
      "column-wise alignment of variable-length loops and are the least reliable here.")
