"""
How many circulating viruses actually carry the surface each design was built against?

Conservation, as this work has used it so far, is a per-position number: how often the
most common residue appears at that position. Two things it cannot tell you, and both
matter more than the thing it can.

The first is whose residue. A design built on BG505 depends on the residues BG505
happens to carry. If BG505 carries a minority variant at a load-bearing position, the
position can be highly conserved and the design still mismatched against most of the
epidemic, because conservation counts the modal residue and the design is not holding
the modal residue. So every hotspot is scored twice here: the frequency of the modal
residue, which is what the paper has been quoting, and the frequency of the residue the
design is actually touching.

The second is coincidence. A binder needs its contacts intact in the same virus, not
each one intact in some virus. Multiplying per-position frequencies would assume the
positions vary independently, which in a recombining virus under linked selection they
do not. The sequences themselves answer it directly: for each of the 2,948 cluster
representatives, ask how much of the design's interface that particular virus presents.

Each contact is weighted by the ddG the alanine scan assigned it, so losing a 4 kcal/mol
anchor counts for more than losing a 1 kcal/mol edge contact. The score per virus is the
fraction of the design's hotspot energy sitting on positions where that virus carries
the same residue as BG505. It is an upper bound on breadth in the sense that a matching
residue is necessary, not sufficient: the surrounding surface may still differ. It is
also not a neutralization prediction.

Two reference sets are scored identically, so the numbers have a scale: the CD4 and 17b
footprints, which are the epitopes of a receptor and of an antibody that really do
engage diverse isolates, and the CCR5 footprint measured in 2.20. For those, every
footprint residue is weighted equally, since no scan exists for them. A size-matched
null is drawn from the scored positions to say what an arbitrary patch of gp120 of the
same size would give.

Run: python strain_coverage.py
"""
import io
import json
import os
import sys
import warnings
from collections import Counter

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_core
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import index_to_one, three_to_index

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = os.path.join(PA, "01_structures")
OUT = os.path.join(PA, "05_tables")

HOT = 1.0            # kcal/mol, the load-bearing threshold used throughout
N_NULL = 2000
SEED = 20261107
FOOTPRINT = 5.0
# V1V2, V3, V4 and V5 vary in length, so a column-wise alignment to HXB2 is weakest
# there and per-position frequencies inside them are the softest numbers here
VAR_LOOPS = [(126, 196), (296, 331), (385, 418), (460, 470)]

DESIGNS = [("D1", "S310"), ("D2", "Chai1_Best"), ("D5", "Latest_Best")]
rng = np.random.default_rng(SEED)
parser = PDBParser(QUIET=True)

M, idx, sizes = conservation_core.matrix()
cons, meta = conservation_core.load()
print("%d cluster representatives, %d scored positions" % (M.shape[0], len(idx)))


def one(res):
    try:
        return index_to_one(three_to_index(res.get_resname()))
    except Exception:
        return None


def column(pos):
    """Residues seen at an HXB2 position across the representatives, gaps dropped."""
    col = M[:, idx[pos]]
    return col, col != 0


def score_set(weights, label, bg505):
    """weights: {hxb2_pos: weight}. bg505: {hxb2_pos: one-letter residue}."""
    pos = [p for p in weights if p in idx and p in bg505]
    if not pos:
        return None
    W = np.array([weights[p] for p in pos], dtype=float)
    tot = W.sum()
    match = np.zeros((M.shape[0], len(pos)), dtype=bool)
    gap = np.zeros_like(match)
    per_pos = []
    for j, p in enumerate(pos):
        col, ok = column(p)
        want = ord(bg505[p])
        match[:, j] = (col == want)
        gap[:, j] = ~ok
        seen = [chr(c) for c in col[ok]]
        c = Counter(seen)
        modal, nmodal = c.most_common(1)[0]
        per_pos.append({"position": p, "bg505_residue": bg505[p],
                        "weight_kcal": round(float(weights[p]), 2),
                        "bg505_frequency": round(float(match[:, j].sum()) / max(ok.sum(), 1), 4),
                        "modal_residue": modal,
                        "modal_frequency": round(nmodal / max(len(seen), 1), 4),
                        "bg505_is_modal": modal == bg505[p],
                        "in_variable_loop": any(lo <= p <= hi for lo, hi in VAR_LOOPS),
                        "target_residue_under_5pct": bool(
                            float(match[:, j].sum()) / max(int(ok.sum()), 1) < 0.05),
                        "gaps": int(gap[:, j].sum())})
    retained = (match * W).sum(1) / tot
    nogap = ~gap.any(1)
    inloop = np.array([q["in_variable_loop"] for q in per_pos])
    # where the shortfall actually sits: weight lost per position, as a share of the
    # whole shortfall, so a diffuse deficit can be told from a localised one
    short = np.array([q["weight_kcal"] * (1 - q["bg505_frequency"]) for q in per_pos])
    for q, v in zip(per_pos, short):
        q["shortfall_kcal"] = round(float(v), 2)
        q["share_of_shortfall_pct"] = round(100 * float(v / short.sum()), 1)             if short.sum() else 0.0
    rank = sorted(per_pos, key=lambda q: -q["shortfall_kcal"])
    out = {"positions_scored": len(pos),
           "total_weight_kcal": round(float(tot), 2),
           "positions_in_variable_loops": int(inloop.sum()),
           "shortfall_kcal": round(float(short.sum()), 2),
           "top5_share_of_shortfall_pct":
               round(sum(q["share_of_shortfall_pct"] for q in rank[:5]), 1),
           "largest_shortfall": ["%s%d %.0f%%" % (q["bg505_residue"], q["position"],
                                                  q["share_of_shortfall_pct"])
                                 for q in rank[:5]],
           "weight_in_variable_loops_pct": round(100 * float(W[inloop].sum() / tot), 1),
           "median_retained_pct": round(100 * float(np.median(retained)), 1),
           "mean_retained_pct": round(100 * float(retained.mean()), 1),
           "strains_with_all_positions_matched_pct":
               round(100 * float((match.all(1)).mean()), 2),
           "strains_retaining_90pct_or_more": round(100 * float((retained >= 0.90).mean()), 1),
           "strains_retaining_75pct_or_more": round(100 * float((retained >= 0.75).mean()), 1),
           "strains_retaining_50pct_or_more": round(100 * float((retained >= 0.50).mean()), 1),
           "gap_free_strains_pct": round(100 * float(nogap.mean()), 1),
           "median_retained_pct_gap_free":
               round(100 * float(np.median(retained[nogap])), 1) if nogap.any() else None,
           "positions_where_bg505_is_not_modal":
               [p["position"] for p in per_pos if not p["bg505_is_modal"]],
           "per_position": per_pos}
    print("%-22s %2d pos | median %5.1f%% retained | all intact in %5.2f%% of strains "
          "| %4.1f%% keep >=75%% | %4.1f%% of weight in a variable loop"
          % (label, len(pos), out["median_retained_pct"],
             out["strains_with_all_positions_matched_pct"],
             out["strains_retaining_75pct_or_more"],
             out["weight_in_variable_loops_pct"]))
    return out


# ---------------------------------------------------------- BG505 residues by position
closedA, openD = {}, {}
for r in next(iter(parser.get_structure("t", os.path.join(SD, "6V8X.pdb"))))["A"]:
    if r.id[0] == " " and one(r):
        closedA[r.id[1]] = one(r)
for r in next(iter(parser.get_structure("o", os.path.join(SD, "6CM3.pdb"))))["D"]:
    if r.id[0] == " " and one(r):
        openD[r.id[1]] = one(r)
# Each contact is scored against the residue in the structure where it was measured:
# the designs were built on 6V8X, the reference footprints are read off 6CM3. The two
# stabilized BG505 constructs are not quite the same molecule, and the disagreements are
# reported rather than averaged away.
bg = dict(openD)
bg.update(closedA)                 # designs: 6V8X preferred
bg_ref = dict(closedA)
bg_ref.update(openD)               # references: 6CM3 preferred
construct_diff = sorted(p for p in set(closedA) & set(openD) if closedA[p] != openD[p])
print("BG505 residues at %d positions; the two constructs differ at %d: %s"
      % (len(bg), len(construct_diff),
         ", ".join("%d (%s in 6V8X, %s in 6CM3)" % (p, closedA[p], openD[p])
                   for p in construct_diff)))

results = {}
for tag, stem in DESIGNS:
    f = os.path.join(OUT, "escape_risk_%s.json" % stem)
    if not os.path.exists(f):
        continue
    d = json.load(io.open(f, encoding="utf-8"))
    w = {c["hxb2_pos"]: c["ddG_kcal_mol"] for c in d["contacts"]
         if c["ddG_kcal_mol"] >= HOT}
    results[tag] = score_set(w, "%s hotspots, ddG wtd" % tag, bg)

# ------------------------------------------------- reference footprints, equal weights
refs = {}
cc = os.path.join(OUT, "ccr5_competition.json")
if os.path.exists(cc):
    d = json.load(io.open(cc, encoding="utf-8"))
    refs["CCR5"] = [int("".join(ch for ch in s if ch.isdigit()))
                    for s in d["ccr5_footprint_on_gp120"]]
    for k in ("CD4", "17b"):
        refs[k] = [int("".join(ch for ch in s if ch.isdigit()))
                   for s in d["partners"][k]["shared"]]
# CD4 and 17b full footprints, measured the same way as everywhere else in this work
op = next(iter(parser.get_structure("o2", os.path.join(SD, "6CM3.pdb"))))
O_res = {r.id[1]: r for r in op["D"] if r.id[0] == " "}


def hv(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def fp_of(chains):
    X = np.vstack([hv(r) for c in chains for r in op[c] if r.id[0] == " " and hv(r).size])
    return [n for n, r in O_res.items() if hv(r).size
            and np.linalg.norm(hv(r)[:, None, :] - X[None, :, :], axis=-1).min()
            < FOOTPRINT]


D_cen = np.vstack([hv(r) for r in op["D"] if r.id[0] == " " and hv(r).size]).mean(0)
cd4_chain = min("GHI", key=lambda c: np.linalg.norm(
    np.vstack([hv(r) for r in op[c] if r.id[0] == " " and hv(r).size]).mean(0) - D_cen))
fab = min((("K", "J"), ("M", "L"), ("O", "N")), key=lambda p: np.linalg.norm(
    np.vstack([hv(r) for c in p for r in op[c] if r.id[0] == " " and hv(r).size]).mean(0)
    - D_cen))
refs["CD4 epitope"] = fp_of([cd4_chain])
refs["17b epitope"] = fp_of(list(fab))
refs.pop("CD4", None)
refs.pop("17b", None)

for name in ("CD4 epitope", "17b epitope", "CCR5"):
    if name in refs:
        results[name] = score_set({p: 1.0 for p in refs[name]},
                                  "%s, equal weights" % name, bg_ref)

# ------------------------------------------------------------------------ null model
# An epitope is a contiguous patch of solvent-exposed surface, so the null has to be
# one too. Scattered positions would be drawn partly from the buried interior, which is
# more conserved than the surface, and the designs would be compared against a standard
# no surface-binding molecule could meet.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
PROBE, NPTS, EXPOSED = 1.4, 128, 20.0      # A, points, A^2 of accessible surface


def _sphere(n=NPTS):
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    theta = np.pi * (1 + 5 ** 0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def residue_sasa(chain):
    """Per-residue solvent accessible area by Shrake-Rupley over the whole chain."""
    from scipy.spatial import cKDTree
    res = [r for r in chain if r.id[0] == " " and one(r)]
    at, owner = [], []
    for i, r in enumerate(res):
        for a in r:
            if a.element != "H":
                at.append(a.coord)
                owner.append((i, VDW.get(a.element, 1.70)))
    X = np.array(at, dtype=float)
    rad = np.array([o[1] for o in owner]) + PROBE
    tree = cKDTree(X)
    S = _sphere()
    area = np.zeros(len(res))
    for i in range(len(X)):
        pts = X[i] + S * rad[i]
        nb = [j for j in tree.query_ball_point(X[i], rad[i] + rad.max()) if j != i]
        if nb:
            d = np.linalg.norm(pts[:, None, :] - X[nb][None, :, :], axis=-1)
            free = (d >= rad[nb][None, :]).all(1)
        else:
            free = np.ones(len(pts), dtype=bool)
        area[owner[i][0]] += 4 * np.pi * rad[i] ** 2 * free.mean()
    return {res[i].id[1]: float(area[i]) for i in range(len(res))}


closed_chain = next(iter(parser.get_structure("s", os.path.join(SD, "6V8X.pdb"))))["A"]
sasa = residue_sasa(closed_chain)
ca = {r.id[1]: r["CA"].coord for r in closed_chain
      if r.id[0] == " " and "CA" in r and one(r)}
surface = sorted(p for p, a in sasa.items()
                 if a >= EXPOSED and p in idx and p in bg and p in ca)
print("\nsurface positions available for the null: %d of %d scored (>=%.0f A^2)"
      % (len(surface), len(idx), EXPOSED))
SX = np.array([ca[p] for p in surface], dtype=float)

nulls = {}
for name in ("D5", "CCR5"):
    if name not in results or not results[name]:
        continue
    n = results[name]["positions_scored"]
    vals = []
    for _ in range(N_NULL):
        c = rng.integers(len(surface))
        order = np.argsort(np.linalg.norm(SX - SX[c], axis=1))[:n]
        pick = [surface[i] for i in order]
        match = np.zeros((M.shape[0], n), dtype=bool)
        for j, p in enumerate(pick):
            col, ok = column(p)
            match[:, j] = (col == ord(bg[p]))
        vals.append(float(np.median(match.sum(1) / n)))
    a = np.array(vals)
    nulls[name] = {"size": n, "draws": N_NULL,
                   "null": "contiguous surface patches of the same size, centred on a "
                           "random exposed residue",
                   "median_of_medians_pct": round(100 * float(np.median(a)), 1),
                   "p5_pct": round(100 * float(np.percentile(a, 5)), 1),
                   "p95_pct": round(100 * float(np.percentile(a, 95)), 1),
                   "percentile_of_observed": round(100 * float(
                       (a <= results[name]["median_retained_pct"] / 100.0).mean()), 1)}
    print("null for %s, %d-residue surface patches: median %.1f%% "
          "(5th-95th %.1f-%.1f); the measured set sits at the %.1fth percentile"
          % (name, n, nulls[name]["median_of_medians_pct"], nulls[name]["p5_pct"],
             nulls[name]["p95_pct"], nulls[name]["percentile_of_observed"]))

# ------------------------------------------- the positions where BG505 is the odd one out
print("\npositions where BG505 does not carry the modal residue:")
for tag in results:
    r = results[tag]
    if not r:
        continue
    odd = [p for p in r["per_position"] if not p["bg505_is_modal"]]
    if odd:
        print("  %-14s %s" % (tag, ", ".join(
            "%s%d (%s %.0f%% vs %s %.0f%%)"
            % (p["bg505_residue"], p["position"], p["bg505_residue"],
               100 * p["bg505_frequency"], p["modal_residue"],
               100 * p["modal_frequency"]) for p in odd)))

with io.open(os.path.join(OUT, "strain_coverage.json"), "w", encoding="utf-8") as fh:
    json.dump({"representatives": int(M.shape[0]),
               "sequences_behind_them": int(sizes.sum()),
               "cluster_identity": meta["cluster_identity"],
               "hotspot_threshold_kcal": HOT,
               "construct_disagreements": construct_diff,
               "residue_source": "designs scored against 6V8X chain A, the model they were built on; reference footprints against 6CM3 chain D, where they were measured",
               "variable_loops": VAR_LOOPS,
               "weighting": "designs weighted by alanine-scan ddG; reference footprints "
                            "weighted equally because no scan exists for them",
               "limits": "a matching residue at a contact position is necessary, not "
                         "sufficient, and this is not a neutralization prediction; "
                         "UniProt's env records are not a balanced sample of the "
                         "epidemic",
               "sets": results, "null_model": nulls}, fh, indent=2)
print("\nwrote strain_coverage.json")
