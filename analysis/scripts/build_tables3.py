"""Tables 16-18: the designs measured in native Env context.

16  accessibility on the closed glycosylated trimer (native_context.json)
17  closed against CD4-bound open state, same isolate (state_dependence.json)
18  rigid-body search for a glycan-free pose (glycan_accommodation.json)
"""
import io
import json
import os

import pandas as pd

TBL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "05_tables")
ORDER = ["CD4", "17b", "D1", "D2", "D2p", "D3", "D5"]
NAME = {"CD4": "CD4 (reference)", "17b": "17b Fab (reference)",
        "D1": "D1 CD4-mimetic control", "D2": "D2 de novo bispecific",
        "D2p": "D2' alternative pose", "D3": "D3 CD4-mimetic cyclic",
        "D5": "D5 de novo, coreceptor site"}


def load(name):
    p = os.path.join(TBL, name)
    if not os.path.exists(p):
        return None
    return json.load(io.open(p, encoding="utf-8"))


def write(df, stem, note=None):
    df.to_csv(os.path.join(TBL, stem + ".csv"), index=False)
    with io.open(os.path.join(TBL, stem + ".md"), "w", encoding="utf-8") as fh:
        fh.write(df.to_markdown(index=False))
        if note:
            fh.write("\n\n" + note + "\n")
    print("wrote", stem)


# ------------------------------------------------------------------- Table 16
nc = load("native_context.json")
if nc:
    rows = []
    for tag in ORDER:
        r = nc["partners"][tag]
        frac = r["ddG_fraction_in_glycan_clash"]
        rows.append({
            "Partner": NAME[tag],
            "Placement": ("native frame" if r["placement"]["frame"] == "native"
                          else "fitted, %.1f A" % r["placement"]["rmsd_A"]),
            "Footprint residues": r["footprint_residues"],
            "Epitope SASA kept (%)": r["epitope_sasa_retained_pct"],
            "Glycan clashes (<3 A)": r["glycan_clashes_3A"],
            "Sequons hit": ", ".join("N%d" % n for n in r["sequons_hit"]) or "none",
            "Neighbour/gp41 clashes (<3 A)": r["neighbour_clashes_3A"],
            "Hotspot energy in glycan contact": ("n/a" if frac is None
                                                 else "%.0f%%" % (100 * frac)),
            "VRC01 Fab clashes (<3 A)": r["vrc01_clashes_3A"],
        })
    write(pd.DataFrame(rows), "table16_native_context",
          "Context is 6V8X (VRC01-bound BG505 F14 SOSIP, 3.0 A). D2, D2', D3 and D5 "
          "carry the deposited protomer coordinates and need no transform; D1, CD4 and "
          "17b are superposed onto their own gp120 footprints and their atom-level "
          "counts inherit that error. Glycans are those resolved in the map, so "
          "occlusion is a lower bound. VRC01 is bound at the CD4bs in this structure, so its "
          "clash column is a check on placement for the CD4bs partners and a test of "
          "whether each design could co-bind a Fab at that site.")

# ------------------------------------------------------------------- Table 17
sd = load("state_dependence.json")
if sd:
    rows = []
    for tag in ("CD4", "17b", "D5", "D2", "D2p", "D3"):
        r = sd["footprints"][tag]
        o = r.get("open_state")
        rows.append({
            "Partner": NAME[tag],
            "Footprint residues": r["footprint_residues"],
            "Resolved in both states": r["resolved_in_both_states"],
            "Local fit RMSD (A)": r.get("local_fit_rmsd_A", "-"),
            "Deviation after protomer fit (A)": r.get("core_fit_deviation_A", "-"),
            "17b clashes in open trimer": o["17b_clashes_3A"] if o else "-",
            "Closest 17b approach (A)": o["17b_min_dist_A"] if o else "-",
            "CD4 clashes in open trimer": o["cd4_clashes_3A"] if o else "-",
            "Glycan clashes in open trimer": o["glycan_clashes_3A"] if o else "-",
            "Neighbour clashes in open trimer": o["neighbour_clashes_3A"] if o else "-",
        })
    write(pd.DataFrame(rows), "table17_state_dependence",
          "Closed state 6V8X chain A, CD4-bound open state 6CM3 chain D, both BG505, "
          "so the comparison carries no isolate or construct difference; whole-protomer "
          "fit between the two states is %.2f A over %d CA. Local fit RMSD is the "
          "footprint superposed on itself and measures whether that surface keeps its "
          "shape through opening. Binders are rigid and are carried onto the open "
          "protomer by their own footprint fit."
          % (sd["whole_protomer_fit_A"], sd["common_CA"]))

# ------------------------------------------------------------------- Table 18
ga = load("glycan_accommodation.json")
if ga:
    rows = []
    for tag in ("D5", "D2", "D2p", "D3"):
        r = ga["partners"].get(tag)
        if not r:
            continue
        m = r["modelled_pose"]
        if r["glycan_free_pose_found"]:
            outcome = "yes"
            move = "%.2f A, %.1f deg" % (r["translation_A"], r["rotation_deg"])
            kept = "%.1f%%" % r["contacts_retained_pct"]
        else:
            outcome = "no"
            move = "none within %.0f A" % r["max_translation_A"]
            kept = "%.1f%% at %d clashes" % (r["contacts_retained_pct"],
                                             r["lowest_glycan_clashes_reached"])
        cur = r.get("clash_budget_curve", {})

        def band(k):
            v = cur.get(k, {}).get("contacts_retained_pct")
            return "not reachable" if v is None else "%.0f%%" % v

        rows.append({
            "Design": NAME[tag],
            "Glycan clashes as modelled": m["glycan_clashes"],
            "gp120 contacts as modelled": m["contacts"],
            "Glycan-free pose found": outcome,
            "Displacement required": move,
            "Interface kept at 0% of original glycan contact": band("0%"),
            "at 10%": band("10%"),
            "at 25%": band("25%"),
            "at 50%": band("50%"),
        })
    write(pd.DataFrame(rows), "table18_glycan_accommodation",
          "Rigid binder and rigid glycans, %d random poses per design within %.0f A and "
          "%.0f degrees of the modelled pose. A pose counts only if it clears every "
          "glycan without pushing into gp120 or a neighbouring protomer. This is an "
          "accessibility search, not an affinity calculation."
          % (ga["search"]["poses"], ga["search"]["max_translation_A"],
             ga["search"]["max_rotation_deg"]))

# ------------------------------------------------------------------- Table 19
ts = load("temsavir_site.json")
if ts:
    rows = []
    for tag in ("D5", "D2", "D2p", "D3"):
        r = ts["designs"].get(tag)
        if not r:
            continue
        rows.append({
            "Design": NAME[tag],
            "Closest approach to temsavir (A)": r["min_dist_to_temsavir_A"],
            "Clashes with temsavir (<3 A)": r["clashes_with_temsavir_3A"],
            "Pocket residues shared": ", ".join(r["pocket_residues_shared"]) or "none",
            "Pocket overlap (%)": r["pocket_overlap_pct"],
        })
    write(pd.DataFrame(rows), "table19_temsavir_site",
          "Temsavir (BMS-626529), the active form of fostemsavir, placed from 5U7O by "
          "superposing its gp120 onto 6V8X chain A (%.2f A over %d CA). The pocket is "
          "the %d gp120 residues within 5 A of the drug."
          % (ts["fit_to_6V8X_A"], ts["common_CA"], len(ts["pocket_residues"])))


# ------------------------------------------------------------------- Table 20
cs = load("counter_screen.json")
if cs:
    rows = []
    for tag, d in cs["designs"].items():
        for tname, t in d["targets"].items():
            rows.append({
                "Binder": d["label"],
                "Target": tname.replace("_", " "),
                "Placements scored": t["n_placements"],
                "Best (kcal/mol)": t["best_kcal"],
                "Median (kcal/mol)": t["median_kcal"],
                "Mean (kcal/mol)": t["mean_kcal"],
                "SD (kcal/mol)": t["sd_kcal"],
                "Designed pose on gp120 (kcal/mol)": d["designed_pose_on_gp120_kcal"],
            })
    write(pd.DataFrame(rows), "table20_counter_screen",
          "Each target screened by an identical procedure: the binder presented face-on "
          "at random surface points, slid to contact, minimised for %d steps and scored "
          "by MM/GBSA. %d placements per target is not a docking search, so each best "
          "score is an upper bound and the bound is loosest for the largest target."
          % (cs["minimisation_steps"], cs["placements_per_target"]))
