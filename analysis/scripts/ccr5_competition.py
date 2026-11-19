"""
Does the design block the actual coreceptor, or only the antibody that stands in for it?

Every competition result in this work so far is measured against 17b, a CD4-induced
antibody used throughout the literature as a surrogate for the coreceptor site. The
claim a reader cares about is the one about CCR5, and the surrogate is not the thing.

6MEO is the gp120-CD4-CCR5 complex at 3.9 A, so the receptor's position on gp120 is
known and can be carried into the frame the designs occupy. Doing that honestly takes
some care, because the two gp120s do not superpose as one rigid body. Over every
residue both resolve the fit is 8 A, and the error is concentrated: the inner domain,
beta20-23 and V5 sit several angstroms apart between the two reconstructions, while the
bridging sheet and outer domain agree to within an angstrom. A global fit would
therefore be the wrong transform to carry a receptor that binds the part which agrees.

So the transfer is local, and it is validated by what it preserves:

  1  6MEO gp120 (chain G) is fitted onto the CD4-bound open protomer 6CM3 chain D over
     the gp120 residues that lie within 10 A of CCR5, which is the surface the receptor
     actually reads. The per-region agreement under the global fit is reported too, so
     the reader can see which parts of the molecule justify a local fit and which do not.
  2  CCR5, with its two sulfotyrosines, is carried in by that transform.
  3  The carry is checked against the thing it has to reproduce: the CCR5-gp120
     interface itself. Contacts and clashes are counted against 6MEO's own gp120 before
     the transfer and against 6CM3 protomer D after it. If those two agree, the
     receptor has landed as a bound receptor and not as an overlapping blob, and that
     is the number that licenses everything after it.

Then the test, with controls that have to come out in opposite directions:

  CD4   binds gp120 at the same time as CCR5 in 6MEO, so it must not clash
  17b   blocks coreceptor engagement functionally, so it must clash
  D5    the de novo design whose whole claim is this site
  D2, D2 prime, D3   the rest of the panel

Each design is carried onto the open protomer by its own footprint fit, exactly as in
2.16, so these numbers inherit that fit error as well as the receptor's 3.9 A map.

CCR5 is a membrane protein. Its transmembrane bundle sits where the membrane is, so
clashes are reported separately for the extracellular N-terminal segment, which carries
the sulfotyrosines that engage the bridging sheet, and for the rest of the receptor.

Run: python ccr5_competition.py
"""
import io
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from scipy.spatial import cKDTree
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import index_to_one, three_to_index

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = os.path.join(PA, "01_structures")
OUT = os.path.join(PA, "05_tables")

FOOTPRINT = 5.0
CLASH, SEVERE, CONTACT = 3.0, 2.5, 4.5
NTERM_LAST = 30          # CCR5 1-30 is the extracellular N-terminal segment
FIT_SHELL = 10.0         # A, gp120 residues near CCR5 define the local superposition
REGIONS = [(31, 125, "inner domain and N terminus"), (126, 196, "V1V2"),
           (197, 250, "bridging sheet and stem"), (251, 295, "outer domain"),
           (296, 331, "V3"), (332, 384, "outer domain, second half"),
           (385, 418, "V4"), (419, 459, "beta20-23"), (460, 491, "V5 and C5")]
GLYCAN = {"NAG", "BMA", "MAN", "A2G"}
KEEP_HET = {"TYS"}       # sulfotyrosine is part of the chain, not a ligand

OPEN_CD4 = ("G", "H", "I")
OPEN_17B = (("K", "J"), ("M", "L"), ("O", "N"))
DESIGNS = [("D5", "Latest_Best_split.pdb"), ("D2", "Chai1_Best_split.pdb"),
           ("D2p", "model_4_chai_split.pdb"), ("D3", "Cyclic_Ring_split.pdb")]

parser = PDBParser(QUIET=True)


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def chain_xyz(chain, glycans=False):
    out = []
    for r in chain:
        ok = (r.id[0] == " " or r.get_resname().strip() in KEEP_HET
              or (glycans and r.get_resname().strip() in GLYCAN))
        if ok:
            out += [a.coord for a in r if a.element != "H"]
    return np.array(out, dtype=float)


def kabsch(mob, ref):
    mc, rc = mob.mean(0), ref.mean(0)
    V, S, Wt = np.linalg.svd((mob - mc).T @ (ref - rc))
    d = np.sign(np.linalg.det(V @ Wt))
    Rm = V @ np.diag([1.0, 1.0, d]) @ Wt
    return Rm, rc - (mc @ Rm)


def fit(mob, ref):
    Rm, tv = kabsch(mob, ref)
    return Rm, tv, float(np.sqrt((((mob @ Rm + tv) - ref) ** 2).sum(-1).mean()))


def one(res):
    try:
        return index_to_one(three_to_index(res.get_resname()))
    except Exception:
        return None


def ca_residues(chain):
    """Residues with a CA, in chain order, with their one-letter codes."""
    out = []
    for r in chain:
        if r.id[0] == " " and "CA" in r and one(r):
            out.append(r)
    return out, "".join(one(r) for r in out)


def correspond(chain_mob, chain_ref):
    """Match residues between two chains by sequence alignment rather than by residue
    number.

    6MEO's gp120 is not numbered in the HXB2 convention the rest of this work uses, and
    the offset is not constant along the chain, so pairing residues by number silently
    compares different residues to each other: over the 362 numbers the two chains share,
    only 28% carry the same amino acid. Aligning the sequences first removes that whole
    class of error, and the identity of the matched pairs is reported so it cannot
    recur unnoticed.
    """
    rm, sm = ca_residues(chain_mob)
    rr, sr = ca_residues(chain_ref)
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"
    a = al.align(sr, sm)[0]
    pairs = []
    for (r0, r1), (m0, m1) in zip(a.aligned[0], a.aligned[1]):
        for k in range(r1 - r0):
            pairs.append((rm[m0 + k], rr[r0 + k]))
    ident = sum(1 for m, r in pairs if one(m) == one(r))
    return pairs, ident


# ------------------------------------------------------------------ structures
op = next(iter(parser.get_structure("o", os.path.join(SD, "6CM3.pdb"))))
me = next(iter(parser.get_structure("m", os.path.join(SD, "6MEO.pdb"))))

O_res = {r.id[1]: r for r in op["D"] if r.id[0] == " " and "CA" in r}
G_res = {r.id[1]: r for r in me["G"] if r.id[0] == " " and "CA" in r}
pairs, ident = correspond(me["G"], op["D"])
print("6MEO gp120 matched to 6CM3 protomer D by sequence: %d residue pairs, "
      "%d identical (%.0f%%)" % (len(pairs), ident, 100.0 * ident / len(pairs)))
MOB = np.array([m["CA"].coord for m, r in pairs], dtype=float)
REF = np.array([r["CA"].coord for m, r in pairs], dtype=float)
ref_num = [r.id[1] for m, r in pairs]

gR, gt, grms = fit(MOB, REF)
dev = np.sqrt((((MOB @ gR + gt) - REF) ** 2).sum(-1))
print("whole-domain fit over all %d pairs: %.2f A" % (len(pairs), grms))
by_region = {}
for lo, hi, lab in REGIONS:
    d = [x for n, x in zip(ref_num, dev) if lo <= n <= hi]
    if d:
        by_region[lab] = {"residues": len(d),
                          "mean_deviation_A": round(float(np.mean(d)), 2)}
        print("   %-28s %3d res  mean deviation %5.2f A" % (lab, len(d), np.mean(d)))

# ------------------------------------------- the surface the receptor reads, fitted local
ccr5_native, nat_map = [], {}
for res in me["B"]:
    if not (res.id[0] == " " or res.get_resname().strip() in KEEP_HET):
        continue
    xyz = heavy(res)
    if xyz.size:
        nat_map[res.id[1]] = (res.get_resname().strip(), xyz)
        ccr5_native.append(xyz)
CCR5_NAT = np.vstack(ccr5_native)
G_NAT = chain_xyz(me["G"])

near = np.array([np.linalg.norm(heavy(m)[:, None, :] - CCR5_NAT[None, :, :],
                                axis=-1).min() < FIT_SHELL for m, r in pairs])
R, t, rms = fit(MOB[near], REF[near])
shell = [n for n, k in zip(ref_num, near) if k]
print("local fit over the %d gp120 residues within %.0f A of CCR5: %.2f A"
      % (int(near.sum()), FIT_SHELL, rms))

D_xyz = chain_xyz(op["D"])
D_cen = D_xyz.mean(0)
cd4_chain = min(OPEN_CD4, key=lambda c: np.linalg.norm(chain_xyz(op[c]).mean(0) - D_cen))
fab = min(OPEN_17B, key=lambda p: np.linalg.norm(
    np.vstack([chain_xyz(op[p[0]]), chain_xyz(op[p[1]])]).mean(0) - D_cen))

# --------------------------------------------------------------- carry the coreceptor
nterm, body, all_ccr5, nmap = [], [], [], {}
for num, (rn, xyz0) in nat_map.items():
    xyz = xyz0 @ R + t
    nmap[num] = (rn, xyz)
    all_ccr5.append(xyz)
    (nterm if num <= NTERM_LAST else body).append(xyz)
CCR5 = np.vstack(all_ccr5)
NT = np.vstack(nterm)
BD = np.vstack(body)
print("CCR5 carried in: %d atoms, %d in the N-terminal segment (1-%d), %d in the body"
      % (len(CCR5), len(NT), NTERM_LAST, len(BD)))

t_ccr5, t_nt, t_bd, t_rec = cKDTree(CCR5), cKDTree(NT), cKDTree(BD), cKDTree(D_xyz)

# the check: the interface the carry has to reproduce
t_gnat = cKDTree(G_NAT)
nat_contacts = int(cKDTree(CCR5_NAT).count_neighbors(t_gnat, CONTACT))
nat_clash = int((t_gnat.query(CCR5_NAT, k=1)[0] < CLASH).sum())
rec_contacts = int(t_ccr5.count_neighbors(t_rec, CONTACT))
rec_clash = int((t_rec.query(CCR5, k=1)[0] < CLASH).sum())
print("CCR5-gp120 interface, native in 6MEO: %d contacts at %.1f A, %d clashes at %.1f A"
      % (nat_contacts, CONTACT, nat_clash, CLASH))
print("                   after the carry: %d contacts at %.1f A, %d clashes at %.1f A"
      % (rec_contacts, CONTACT, rec_clash, CLASH))

ccr5_fp = sorted(n for n, r in O_res.items() if heavy(r).size
                 and np.linalg.norm(heavy(r)[:, None, :] - CCR5[None, :, :],
                                    axis=-1).min() < FOOTPRINT)
reach = sorted(n for n, (rn, xyz) in nmap.items()
               if np.linalg.norm(xyz[:, None, :] - D_xyz[None, :, :], axis=-1).min()
               < FOOTPRINT)
tys = [n for n in reach if nmap[n][0] == "TYS"]
print("CCR5 footprint on gp120: %d residues | %d CCR5 residues reach gp120, "
      "sulfotyrosines among them: %s"
      % (len(ccr5_fp), len(reach), ", ".join("Tys%d" % n for n in tys) or "none"))


# ------------------------------------------------------------------ the comparison
def vs(tree, P, cut=CLASH):
    d = tree.query(P, k=1)[0]
    return int((d < cut).sum()), round(float(d.min()), 2)


def footprint_of(xyz):
    return sorted(n for n, r in O_res.items() if heavy(r).size
                  and np.linalg.norm(heavy(r)[:, None, :] - xyz[None, :, :],
                                     axis=-1).min() < FOOTPRINT)


def row_for(P, fp):
    ov = sorted(set(fp) & set(ccr5_fp))
    nc, nmin = vs(t_nt, P)
    bc, bmin = vs(t_bd, P)
    return {"footprint_residues": len(fp),
            "ccr5_clashes_3A": vs(t_ccr5, P)[0],
            "ccr5_severe_2.5A": vs(t_ccr5, P, SEVERE)[0],
            "ccr5_min_dist_A": vs(t_ccr5, P)[1],
            "contacts_with_ccr5_4.5A": int(cKDTree(P).count_neighbors(t_ccr5, CONTACT)),
            "nterm_clashes_3A": nc, "nterm_min_dist_A": nmin,
            "body_clashes_3A": bc, "body_min_dist_A": bmin,
            "shared_footprint_residues": len(ov),
            "pct_of_ccr5_footprint_covered": round(100.0 * len(ov) / len(ccr5_fp), 1),
            "shared": ["%s%d" % (O_res[n].get_resname().title(), n) for n in ov]}


rows = {}
CD4_O = chain_xyz(op[cd4_chain])
FAB_O = np.vstack([chain_xyz(op[fab[0]]), chain_xyz(op[fab[1]])])
rows["CD4"] = row_for(CD4_O, footprint_of(CD4_O))
rows["CD4"]["note"] = ("negative control: CD4 and CCR5 are bound to the same gp120 in "
                       "6MEO, so this must not clash")
rows["17b"] = row_for(FAB_O, footprint_of(FAB_O))
rows["17b"]["note"] = ("positive control: 17b blocks coreceptor engagement, so this "
                       "must clash")

for tag, fname in DESIGNS:
    st = next(iter(parser.get_structure("d", os.path.join(SD, fname))))
    B = np.array([a.coord for r in st["B"] if r.id[0] == " "
                  for a in r if a.element != "H"], dtype=float)
    C_res = {r.id[1]: r for r in st["A"] if r.id[0] == " " and "CA" in r}
    fp_closed = [n for n, r in C_res.items() if heavy(r).size
                 and np.linalg.norm(heavy(r)[:, None, :] - B[None, :, :],
                                    axis=-1).min() < FOOTPRINT]
    use = [n for n in sorted(fp_closed) if n in O_res]
    lR, lt, lrms = fit(np.array([C_res[n]["CA"].coord for n in use], dtype=float),
                       np.array([O_res[n]["CA"].coord for n in use], dtype=float))
    P = B @ lR + lt
    r = row_for(P, footprint_of(P))
    r["placed_by"] = "own footprint fit, %d CA at %.2f A" % (len(use), lrms)
    rows[tag] = r

hdr = ("%-5s %4s %7s %7s %8s %8s %8s %7s %7s"
       % ("tag", "fp", "clash", "severe", "minDist", "contact", "Nterm", "body", "%CCR5"))
print("\n" + hdr)
print("-" * len(hdr))
for tag in ("CD4", "17b", "D5", "D2", "D2p", "D3"):
    r = rows[tag]
    print("%-5s %4d %7d %7d %8.2f %8d %8d %7d %7.1f"
          % (tag, r["footprint_residues"], r["ccr5_clashes_3A"], r["ccr5_severe_2.5A"],
             r["ccr5_min_dist_A"], r["contacts_with_ccr5_4.5A"], r["nterm_clashes_3A"],
             r["body_clashes_3A"], r["pct_of_ccr5_footprint_covered"]))

# ------------------------------- is escape at the design's hotspots costly to the virus?
hot = {}
for tag, stem in (("D5", "Latest_Best"), ("D2", "Chai1_Best"), ("D1", "S310")):
    f = os.path.join(OUT, "escape_risk_%s.json" % stem)
    if not os.path.exists(f):
        continue
    d = json.load(io.open(f, encoding="utf-8"))
    cs = [c for c in d["contacts"] if c["ddG_kcal_mol"] >= 1.0]
    on = [c for c in cs if c["hxb2_pos"] in set(ccr5_fp)]
    tot = sum(c["ddG_kcal_mol"] for c in cs) or 1.0
    hot[tag] = {"load_bearing_contacts": len(cs),
                "also_contacted_by_ccr5": len(on),
                "pct_of_hotspot_energy_on_ccr5_contacts":
                    round(100.0 * sum(c["ddG_kcal_mol"] for c in on) / tot, 1),
                "positions": ["%s%d" % (c["resname"].title(), c["hxb2_pos"])
                              for c in on]}
    print("%s: %d of %d load-bearing contacts are also CCR5 contacts, %.1f%% of its "
          "hotspot energy" % (tag, len(on), len(cs),
                              hot[tag]["pct_of_hotspot_energy_on_ccr5_contacts"]))

with io.open(os.path.join(OUT, "ccr5_competition.json"), "w", encoding="utf-8") as fh:
    json.dump({"coreceptor_source": "6MEO, gp120-CD4-CCR5, cryo-EM 3.90 A",
               "open_trimer": "6CM3 chain D, BG505 SOSIP with sCD4, 17b and 8ANC195",
               "global_fit_A": round(grms, 2), "matched_pairs": len(pairs), "identical_in_matched_pairs": ident,
               "global_fit_by_region": by_region,
               "local_fit_shell_A": FIT_SHELL, "local_fit_residues": int(near.sum()),
               "local_fit_A": round(rms, 2),
               "carry_validation": {
                   "note": "the CCR5-gp120 interface before and after the transfer; the "
                           "carry is only usable if it reproduces the interface it came "
                           "from",
                   "native_contacts_4.5A": nat_contacts,
                   "carried_contacts_4.5A": rec_contacts,
                   "native_clashes_3A": nat_clash,
                   "carried_clashes_3A": rec_clash},
               "ccr5_atoms": int(len(CCR5)), "nterm_last_residue": NTERM_LAST,
               "ccr5_footprint_on_gp120": ["%s%d" % (O_res[n].get_resname().title(), n)
                                           for n in ccr5_fp],
               "ccr5_residues_reaching_gp120": len(reach),
               "sulfotyrosines_in_contact": ["Tys%d" % n for n in tys],
               "caveat": "the coreceptor is carried by a cross-isolate gp120 fit and the "
                         "designs by their own footprint fits, so every number inherits "
                         "both; CCR5 and the designs are rigid",
               "partners": rows, "hotspot_overlap": hot}, fh, indent=2)
print("\nwrote ccr5_competition.json")
