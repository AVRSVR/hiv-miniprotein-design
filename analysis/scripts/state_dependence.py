"""
Does the epitope survive CD4-induced opening, and can the binder still reach it?

17b binds only after CD4 has opened the trimer, because its epitope forms on opening.
That is the central liability of CD4i-directed agents, and D5 blocks the same site, so
the question applies to D5 directly: is the surface it binds present in both states?

Both states of the same isolate are available, which removes the isolate and construct
differences that make a comparison against the 1GC1 HXB2 core uninterpretable:

  closed  6V8X chain A, VRC01-bound BG505 F14 SOSIP, the model every design was built on
  open    6CM3 chain D, BG505 SOSIP with soluble CD4, 17b and 8ANC195 bound

Per footprint:
  local fit RMSD   superposing the two states over the footprint residues alone, which
                   measures whether that surface keeps its own shape through the
                   transition, independent of any rigid-body movement
  core deviation   where the footprint lands once the protomers are aligned globally

Each binder is then carried onto the open protomer by its own footprint fit and measured
against the open trimer: its glycans, its neighbouring protomers and gp41, and the CD4
and 17b molecules actually bound to that protomer. This re-tests the paper's central
competition claim in the native open trimer rather than by superposition onto the
1GC1 core.

Run: python state_dependence.py
"""
import io
import json
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = os.path.join(PA, "01_structures")
OUT = os.path.join(PA, "05_tables")

FOOTPRINT = 5.0
CLASH, SEVERE, CONTACT = 3.0, 2.5, 4.5
GLYCAN = {"NAG", "BMA", "MAN"}

# 6CM3 chain roles
OPEN_GP120 = ("D", "E", "F")
OPEN_GP41 = ("A", "B", "C")
OPEN_CD4 = ("G", "H", "I")
OPEN_17B = (("K", "J"), ("M", "L"), ("O", "N"))      # heavy, light
OPEN_8ANC = ("P", "Q", "R", "S", "T", "U")

DESIGNS = [("D5", "Latest_Best_split.pdb"), ("D2", "Chai1_Best_split.pdb"),
           ("D2p", "model_4_chai_split.pdb"), ("D3", "Cyclic_Ring_split.pdb")]

parser = PDBParser(QUIET=True)


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def heavy_chain(chain, hetero=False):
    return np.array([a.coord for r in chain
                     if (r.id[0] == " " or (hetero and r.get_resname().strip() in GLYCAN))
                     for a in r if a.element != "H"], dtype=float)


def kabsch(mob, ref):
    mc, rc = mob.mean(0), ref.mean(0)
    V, S, Wt = np.linalg.svd((mob - mc).T @ (ref - rc))
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, rc - (mc @ R)


def rmsd_after(mob, ref):
    R, t = kabsch(mob, ref)
    return R, t, float(np.sqrt((((mob @ R + t) - ref) ** 2).sum(-1).mean()))


closed = next(iter(parser.get_structure("c", os.path.join(SD, "6V8X.pdb"))))
op = next(iter(parser.get_structure("o", os.path.join(SD, "6CM3.pdb"))))

C_res = {r.id[1]: r for r in closed["A"] if r.id[0] == " " and "CA" in r}
O_res = {r.id[1]: r for r in op["D"] if r.id[0] == " " and "CA" in r}
common = sorted(set(C_res) & set(O_res))
core_R, core_t, core_rms = rmsd_after(
    np.array([C_res[n]["CA"].coord for n in common], dtype=float),
    np.array([O_res[n]["CA"].coord for n in common], dtype=float))
print("BG505 closed (6V8X A) vs CD4-bound open (6CM3 D): %d common CA, "
      "whole-protomer fit %.2f A" % (len(common), core_rms))

# ------------------------------------------------------------- open trimer context
gly_o, neigh_o = [], []
for ch in op:
    for res in ch:
        if res.get_resname().strip() in GLYCAN:
            gly_o += [a.coord for a in res if a.element != "H"]
        elif res.id[0] == " " and (ch.id in OPEN_GP41
                                   or (ch.id in OPEN_GP120 and ch.id != "D")):
            neigh_o += [a.coord for a in res if a.element != "H"]
GLY_O = np.array(gly_o, dtype=float)
NEIGH_O = np.array(neigh_o, dtype=float)
D_xyz = heavy_chain(op["D"])

# the CD4 and 17b copies bound to protomer D
D_cen = D_xyz.mean(0)
cd4_chain = min(OPEN_CD4, key=lambda c: np.linalg.norm(heavy_chain(op[c]).mean(0) - D_cen))
fab = min(OPEN_17B, key=lambda p: np.linalg.norm(
    np.vstack([heavy_chain(op[p[0]]), heavy_chain(op[p[1]])]).mean(0) - D_cen))
CD4_O = heavy_chain(op[cd4_chain])
FAB_O = np.vstack([heavy_chain(op[fab[0]]), heavy_chain(op[fab[1]])])
print("protomer D is engaged by CD4 chain %s and 17b chains %s/%s; "
      "%d glycan atoms, %d neighbouring protomer and gp41 atoms"
      % (cd4_chain, fab[0], fab[1], len(GLY_O), len(NEIGH_O)))

t_gly, t_nb, t_rec, t_cd4, t_fab = (cKDTree(GLY_O), cKDTree(NEIGH_O), cKDTree(D_xyz),
                                    cKDTree(CD4_O), cKDTree(FAB_O))


def vs(tree, P, cut=CLASH):
    d = tree.query(P, k=1)[0]
    return int((d < cut).sum()), round(float(d.min()), 2)


# ----------------------------------------------------------------- per footprint
def design_pose(fname):
    st = next(iter(parser.get_structure("d", os.path.join(SD, fname))))
    B = np.array([a.coord for r in st["B"] if r.id[0] == " "
                  for a in r if a.element != "H"], dtype=float)
    g = st["A"]
    fp = [r.id[1] for r in g if r.id[0] == " " and heavy(r).size
          and np.linalg.norm(heavy(r)[:, None, :] - B[None, :, :], axis=-1).min()
          < FOOTPRINT]
    return B, fp


def ligand_footprint_closed(xyz):
    return [n for n, r in C_res.items() if heavy(r).size
            and np.linalg.norm(heavy(r)[:, None, :] - xyz[None, :, :], axis=-1).min()
            < FOOTPRINT]


cases = {}
for tag, fname in DESIGNS:
    B, fp = design_pose(fname)
    cases[tag] = (B, fp)
# the two experimental partners, footprints read off the open structure where they are
# actually bound, then compared between states
cases["CD4"] = (None, [n for n, r in O_res.items() if heavy(r).size
                       and np.linalg.norm(heavy(r)[:, None, :] - CD4_O[None, :, :],
                                          axis=-1).min() < FOOTPRINT])
cases["17b"] = (None, [n for n, r in O_res.items() if heavy(r).size
                       and np.linalg.norm(heavy(r)[:, None, :] - FAB_O[None, :, :],
                                          axis=-1).min() < FOOTPRINT])

rows = {}
for tag, (B, fp) in cases.items():
    use = [n for n in sorted(fp) if n in C_res and n in O_res]
    row = {"footprint_residues": len(fp), "resolved_in_both_states": len(use),
           "unresolved_in_one_state": len(fp) - len(use)}
    if len(use) >= 4:
        Cc = np.array([C_res[n]["CA"].coord for n in use], dtype=float)
        Oc = np.array([O_res[n]["CA"].coord for n in use], dtype=float)
        row["core_fit_deviation_A"] = round(
            float(np.sqrt((((Cc @ core_R + core_t) - Oc) ** 2).sum(-1).mean())), 2)
        lR, lt, lrms = rmsd_after(Cc, Oc)
        row["local_fit_rmsd_A"] = round(lrms, 2)
    if B is not None and len(use) >= 4:
        P = B @ lR + lt                      # carry the binder on its own footprint fit
        row["open_state"] = {
            "placed_by": "footprint fit, %d CA at %.2f A" % (len(use), lrms),
            "receptor_clashes_2.5A": vs(t_rec, P, SEVERE)[0],
            "glycan_clashes_3A": vs(t_gly, P)[0],
            "glycan_min_dist_A": vs(t_gly, P)[1],
            "neighbour_clashes_3A": vs(t_nb, P)[0],
            "neighbour_min_dist_A": vs(t_nb, P)[1],
            "cd4_clashes_3A": vs(t_cd4, P)[0],
            "cd4_min_dist_A": vs(t_cd4, P)[1],
            "17b_clashes_3A": vs(t_fab, P)[0],
            "17b_min_dist_A": vs(t_fab, P)[1],
            "contacts_with_receptor_4.5A": int(cKDTree(P).count_neighbors(t_rec, CONTACT)),
        }
    rows[tag] = row

hdr = ("%-5s %4s %5s %9s %9s | %7s %7s %7s %7s %7s"
       % ("tag", "fp", "both", "coreDev", "localRMS", "17bCl", "17bMin", "CD4Cl",
          "glyCl", "nbrCl"))
print("\n" + hdr)
print("-" * len(hdr))
for tag in ("CD4", "17b", "D5", "D2", "D2p", "D3"):
    r = rows[tag]
    o = r.get("open_state")
    print("%-5s %4d %5d %9s %9s | %7s %7s %7s %7s %7s"
          % (tag, r["footprint_residues"], r["resolved_in_both_states"],
             r.get("core_fit_deviation_A", "-"), r.get("local_fit_rmsd_A", "-"),
             o["17b_clashes_3A"] if o else "-", o["17b_min_dist_A"] if o else "-",
             o["cd4_clashes_3A"] if o else "-", o["glycan_clashes_3A"] if o else "-",
             o["neighbour_clashes_3A"] if o else "-"))

with io.open(os.path.join(OUT, "state_dependence.json"), "w", encoding="utf-8") as fh:
    json.dump({"closed_state": "6V8X chain A, VRC01-bound BG505 F14 SOSIP",
               "open_state": "6CM3 chain D, BG505 SOSIP with sCD4, 17b and 8ANC195",
               "same_isolate": True,
               "whole_protomer_fit_A": round(core_rms, 2), "common_CA": len(common),
               "cd4_chain": cd4_chain, "fab_chains": list(fab),
               "caveat": "binders are rigid and are carried onto the open protomer by "
                         "their own footprint fit, so their open-state numbers inherit "
                         "that fit error; glycans are those resolved in each map",
               "footprints": rows}, fh, indent=2)
print("\nwrote state_dependence.json")
