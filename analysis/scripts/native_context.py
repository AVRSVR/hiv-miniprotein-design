"""
Can these binders reach their epitopes on the native Env trimer?

Every interface in this study was computed on a single gp120 protomer with the
heteroatoms stripped. On the virion gp120 is one subunit of a glycosylated trimer, so a
footprint that looks open on an isolated protomer may be covered by the neighbouring
protomer, by gp41, or by a glycan.

The design target was extracted from 6V8X (VRC01-bound BG505 F14 SOSIP trimer, 3.0 A)
and four of the five complexes still carry that protomer at its deposited coordinates
(chain A CA RMSD 0.00 A against 6V8X chain A). The rest of the deposited structure can
therefore be dropped back in with no superposition at all, which restores the native
closed prefusion trimer, its gp41 subunits, its modelled glycans and the VRC01 Fab
bound at the CD4bs.

The remaining complex, D1, carries a different gp120 model: it fits the CD4-bound cores
of 1RZJ and 1GC1 at 1.75 and 2.33 A but 6V8X chain A at 16.5 A over 434 CA, so it is in
the CD4-bound rather than the closed conformation. D1, CD4 and the 17b Fab are therefore
placed by Kabsch superposition over their own gp120 footprints, and their placement
error is reported alongside their clash counts.

Measured per partner:
  clashes           partner heavy atoms within 2.5 / 3.0 / 4.5 A of trimer context
  which glycans     the sequon asparagines whose glycans are hit
  which residues    the partner's own clashing residues, with their alanine-scan ddG,
                    which says whether glycan contact falls on load-bearing residues
  accessibility     per-residue SASA of the gp120 footprint in the isolated protomer
                    against the same residues in the assembled glycosylated trimer

Run: python native_context.py
"""
import io
import json
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = os.path.join(PA, "01_structures")
OUT = os.path.join(PA, "05_tables")

TRIMER = os.path.join(SD, "6V8X.pdb")
GLYCAN = {"NAG", "BMA", "MAN"}
GP120 = ("A", "E", "I")          # A is the design target protomer
GP41 = ("B", "F", "J")
VRC01 = ("C", "D", "G", "H", "K", "L")

DESIGNS = [("D1", "S310_split.pdb", "alascan_S310.json"),
           ("D2", "Chai1_Best_split.pdb", "alascan_Chai1_Best.json"),
           ("D2p", "model_4_chai_split.pdb", None),
           ("D3", "Cyclic_Ring_split.pdb", "alascan_Cyclic_Ring.json"),
           ("D5", "Latest_Best_split.pdb", "alascan_Latest_Best.json")]

FOOTPRINT = 5.0        # A, gp120 residue counts as footprint within this of partner
SEVERE, CLASH, NEAR = 2.5, 3.0, 4.5
NATIVE_TOL = 0.5       # A, chain A RMSD below this means the complex needs no transform

parser = PDBParser(QUIET=True)


def heavy_coords(res_iter):
    return np.array([a.coord for r in res_iter if r.id[0] == " "
                     for a in r if a.element != "H"], dtype=float)


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def kabsch(mob, ref):
    mc, rc = mob.mean(0), ref.mean(0)
    V, S, Wt = np.linalg.svd((mob - mc).T @ (ref - rc))
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, rc - (mc @ R)


def contacts(A, B, chunk=3000):
    """Minimum distance, clash counts, and the indices of A that clash."""
    best = np.inf
    n_sev = n_cl = n_near = 0
    hits = set()
    for i in range(0, len(A), chunk):
        D = np.linalg.norm(A[i:i + chunk, None, :] - B[None, :, :], axis=-1)
        best = min(best, float(D.min()))
        n_sev += int((D < SEVERE).sum())
        n_cl += int((D < CLASH).sum())
        n_near += int((D < NEAR).sum())
        hits |= {i + int(k) for k in np.where((D < CLASH).any(axis=1))[0]}
    return best, n_sev, n_cl, n_near, hits


# ---------------------------------------------------------------- trimer context
tri = next(iter(parser.get_structure("t", TRIMER)))

glyc_atoms, glyc_res, prot_other, gp41_atoms, vrc01_atoms = [], [], [], [], []
for ch in tri:
    for res in ch:
        if res.get_resname().strip() in GLYCAN:
            at = [a for a in res if a.element != "H"]
            glyc_atoms += at
            glyc_res.append((ch.id, res.id[1], at))
            continue
        if res.id[0] != " ":
            continue
        if ch.id in GP41:
            gp41_atoms += [a for a in res if a.element != "H"]
        elif ch.id in VRC01:
            vrc01_atoms += [a for a in res if a.element != "H"]
        elif ch.id in GP120 and ch.id != "A":
            prot_other += [a for a in res if a.element != "H"]

GLY = np.array([a.coord for a in glyc_atoms], dtype=float)
NEIGH = np.array([a.coord for a in prot_other + gp41_atoms], dtype=float)
VRC = np.array([a.coord for a in vrc01_atoms], dtype=float)
print("context: %d glycan atoms in %d sugar residues, %d neighbouring protomer and "
      "gp41 atoms, %d VRC01 atoms" % (len(GLY), len(glyc_res), len(NEIGH), len(VRC)))

prot_A = [a for r in tri["A"] if r.id[0] == " " for a in r if a.element != "H"]
A_res = {r.id[1]: r for r in tri["A"] if r.id[0] == " "}
A_ca = {n: r["CA"].coord for n, r in A_res.items() if "CA" in r}

# each glycan residue's nearest chain A asparagine, so clashes can be named by sequon
asn = {n: r for n, r in A_res.items() if r.get_resname() == "ASN"}
asn_xyz = {n: heavy(r) for n, r in asn.items()}
glyc_sequon = []
for cid, rid, at in glyc_res:
    G = np.array([a.coord for a in at], dtype=float)
    best, bn = np.inf, None
    for n, X in asn_xyz.items():
        d = float(np.linalg.norm(G[:, None, :] - X[None, :, :], axis=-1).min())
        if d < best:
            best, bn = d, n
    glyc_sequon.append(bn if best < 2.0 else None)


# ------------------------------------------------------------ partners to measure
def load_partner(fname, binder_chain="B", gp120_chain="A"):
    """Binder heavy atoms, plus the transform that puts them in the 6V8X frame."""
    st = next(iter(parser.get_structure("d", os.path.join(SD, fname))))
    bres = [r for r in st[binder_chain] if r.id[0] == " "]
    B = heavy_coords(bres)
    own = {r.id[1]: r["CA"].coord for r in st[gp120_chain]
           if r.id[0] == " " and "CA" in r}
    com = sorted(set(own) & set(A_ca))
    M = np.array([own[i] for i in com], dtype=float)
    R0 = np.array([A_ca[i] for i in com], dtype=float)
    raw = float(np.sqrt(((M - R0) ** 2).sum(-1).mean()))
    if raw < NATIVE_TOL:
        return B, bres, {"frame": "native", "rmsd_A": round(raw, 2)}
    # local fit over the gp120 residues this binder actually touches
    near = [i for i in com
            if np.linalg.norm(heavy(st[gp120_chain][i])[:, None, :]
                              - B[None, :, :], axis=-1).min() < 8.0]
    if len(near) < 8:
        near = com
    M2 = np.array([own[i] for i in near], dtype=float)
    R2 = np.array([A_ca[i] for i in near], dtype=float)
    R, t = kabsch(M2, R2)
    rms = float(np.sqrt((((M2 @ R + t) - R2) ** 2).sum(-1).mean()))
    return B @ R + t, bres, {"frame": "fitted", "epitope_CA": len(near),
                             "rmsd_A": round(rms, 2)}


def load_reference(ligand_chains):
    """Place a 1GC1 ligand into the 6V8X frame by fitting its own gp120 epitope."""
    st = next(iter(parser.get_structure("r", os.path.join(SD, "1GC1.pdb"))))
    g = st["G"]
    L = heavy_coords([r for c in ligand_chains for r in st[c]])
    lres = [r for c in ligand_chains for r in st[c] if r.id[0] == " "]
    epi = [r.id[1] for r in g if r.id[0] == " " and "CA" in r
           and heavy(r).size
           and np.linalg.norm(heavy(r)[:, None, :] - L[None, :, :], axis=-1).min()
           < FOOTPRINT]
    pairs = [n for n in epi if n in A_ca and "CA" in g[n]]
    M = np.array([g[n]["CA"].coord for n in pairs], dtype=float)
    R0 = np.array([A_ca[n] for n in pairs], dtype=float)
    R, t = kabsch(M, R0)
    rms = float(np.sqrt((((M @ R + t) - R0) ** 2).sum(-1).mean()))
    return L @ R + t, lres, {"frame": "fitted", "epitope_CA": len(pairs),
                             "rmsd_A": round(rms, 2)}


partners = {}
for tag, fname, ala in DESIGNS:
    B, bres, place = load_partner(fname)
    ddg = {}
    if ala and os.path.exists(os.path.join(OUT, ala)):
        ddg = {r["resid"]: r["ddG_kcal_mol"]
               for r in json.load(io.open(os.path.join(OUT, ala), encoding="utf-8"))["residues"]}
    partners[tag] = (B, bres, place, ddg)
    print("%-4s %s frame, placement RMSD %.2f A" % (tag, place["frame"], place["rmsd_A"]))
for tag, chains in (("CD4", ("C",)), ("17b", ("H", "L"))):
    B, bres, place = load_reference(chains)
    partners[tag] = (B, bres, place, {})
    print("%-4s placed on %d epitope CA, fit RMSD %.2f A"
          % (tag, place["epitope_CA"], place["rmsd_A"]))


# --------------------------------------------------- epitope accessibility by SASA
# Shrake-Rupley over a flat atom list. Biopython's implementation only accepts whole
# entities, and the assembled system here mixes chains and sugars from several of them.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
PROBE = 1.4
NPTS = 128


def _sphere(n=NPTS):
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    theta = np.pi * (1 + 5 ** 0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


SPHERE = _sphere()


def sasa_of(atoms):
    """Per-atom SASA for a flat atom list; returns dict(id(atom) -> area)."""
    from scipy.spatial import cKDTree
    X = np.array([a.coord for a in atoms], dtype=float)
    R = np.array([VDW.get(a.element, 1.70) for a in atoms], dtype=float) + PROBE
    tree = cKDTree(X)
    rmax = R.max()
    out = {}
    for i, a in enumerate(atoms):
        pts = X[i] + R[i] * SPHERE
        nb = [j for j in tree.query_ball_point(X[i], R[i] + rmax) if j != i]
        if nb:
            D = np.linalg.norm(pts[:, None, :] - X[nb][None, :, :], axis=-1)
            free = np.all(D >= R[nb][None, :], axis=1)
        else:
            free = np.ones(len(pts), bool)
        out[id(a)] = 4.0 * np.pi * R[i] ** 2 * free.mean()
    return out


iso = sasa_of(prot_A)
full = sasa_of(prot_A + prot_other + gp41_atoms + glyc_atoms)

rows = {}
for tag, (P, bres, place, ddg) in partners.items():
    fp = sorted(n for n, r in A_res.items() if heavy(r).size
                and np.linalg.norm(heavy(r)[:, None, :] - P[None, :, :], axis=-1).min()
                < FOOTPRINT)
    per_res = []
    for n in fp:
        at = [a for a in A_res[n] if a.element != "H"]
        per_res.append((n, A_res[n].get_resname(),
                        sum(iso[id(a)] for a in at), sum(full[id(a)] for a in at)))
    tot_iso = sum(p[2] for p in per_res) or 1.0
    tot_full = sum(p[3] for p in per_res)
    occluded = [p for p in per_res if p[2] > 1.0 and p[3] / p[2] < 0.25]

    g_min, g_sev, g_cl, g_near, g_hits = contacts(P, GLY)
    n_min, n_sev, n_cl, _, n_hits = contacts(P, NEIGH)
    v_min, _, v_cl, _, _ = contacts(P, VRC)

    # which sequons are hit, and which of the partner's own residues do the hitting
    hit_sequons = set()
    for k, (cid, rid, at) in enumerate(glyc_res):
        G = np.array([a.coord for a in at], dtype=float)
        if np.linalg.norm(P[:, None, :] - G[None, :, :], axis=-1).min() < CLASH:
            if glyc_sequon[k] is not None:
                hit_sequons.add(glyc_sequon[k])

    # map clashing partner atoms back to partner residues
    idx, res_of_atom = 0, []
    for r in bres:
        for a in r:
            if a.element != "H":
                res_of_atom.append(r)
                idx += 1
    clash_res = sorted({(res_of_atom[i].id[1], res_of_atom[i].get_resname())
                        for i in g_hits if i < len(res_of_atom)})
    clash_ddg = round(sum(max(ddg.get(n, 0.0), 0.0) for n, _ in clash_res), 1)
    tot_ddg = round(sum(v for v in ddg.values() if v > 0), 1) if ddg else None

    rows[tag] = {
        "placement": place,
        "footprint_residues": len(fp),
        "epitope_sasa_isolated_A2": round(tot_iso, 1),
        "epitope_sasa_in_trimer_A2": round(tot_full, 1),
        "epitope_sasa_retained_pct": round(100.0 * tot_full / tot_iso, 1),
        "footprint_residues_occluded": len(occluded),
        "occluded_residues": [f"{p[1]}{p[0]}" for p in occluded],
        "glycan_min_dist_A": round(g_min, 2),
        "glycan_clashes_2.5A": g_sev,
        "glycan_clashes_3A": g_cl,
        "glycan_contacts_4.5A": g_near,
        "sequons_hit": sorted(hit_sequons),
        "own_residues_in_glycan_clash": [f"{x[1]}{x[0]}" for x in clash_res],
        "ddG_of_clashing_residues_kcal": clash_ddg,
        "ddG_total_hotspots_kcal": tot_ddg,
        "ddG_fraction_in_glycan_clash": (round(clash_ddg / tot_ddg, 3)
                                         if tot_ddg else None),
        "neighbour_min_dist_A": round(n_min, 2),
        "neighbour_clashes_2.5A": n_sev,
        "neighbour_clashes_3A": n_cl,
        "vrc01_min_dist_A": round(v_min, 2),
        "vrc01_clashes_3A": v_cl,
    }

hdr = ("%-5s %-7s %4s %7s %6s %7s %6s %7s %6s %7s"
       % ("tag", "frame", "fpr", "keep%", "occl", "glyMin", "glyCl", "nbrMin",
          "nbrCl", "ddGhit"))
print("\n" + hdr)
print("-" * len(hdr))
for tag in ("CD4", "17b", "D1", "D2", "D2p", "D3", "D5"):
    r = rows[tag]
    print("%-5s %-7s %4d %7.1f %6d %7.2f %6d %7.2f %6d %7s"
          % (tag, r["placement"]["frame"], r["footprint_residues"],
             r["epitope_sasa_retained_pct"], r["footprint_residues_occluded"],
             r["glycan_min_dist_A"], r["glycan_clashes_3A"],
             r["neighbour_min_dist_A"], r["neighbour_clashes_3A"],
             r["ddG_fraction_in_glycan_clash"]))
print("\nsequons hit:")
for tag in ("CD4", "17b", "D1", "D2", "D2p", "D3", "D5"):
    print("  %-5s %s" % (tag, rows[tag]["sequons_hit"] or "none"))

meta = {"trimer": "6V8X, VRC01-bound BG505 F14 SOSIP trimer, 3.0 A",
        "frame": "D2/D2p/D3/D5 are in the deposited 6V8X frame with no transform; "
                 "D1, CD4 and 17b are superposed onto their own gp120 footprints",
        "glycan_atoms": int(len(GLY)), "glycan_residues": len(glyc_res),
        "neighbour_atoms": int(len(NEIGH)),
        "caveat": "glycans are those resolved in the cryo-EM map, which model the core "
                  "and under-represent the processed glycan, so occlusion is a lower bound",
        "partners": rows}
with io.open(os.path.join(OUT, "native_context.json"), "w", encoding="utf-8") as fh:
    json.dump(meta, fh, indent=2)
print("\nwrote native_context.json")
