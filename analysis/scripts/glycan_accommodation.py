"""
Is the glycan conflict fatal, or can the binder shift out of it?

native_context.py measures each partner where it was modelled and finds that most of
them clash with the glycans resolved in 6V8X. A clash in a fixed pose does not settle
the question: a binder that can shed its glycan contacts with a small rigid-body
adjustment while keeping most of its protein contacts is still a candidate, whereas one
whose protein contacts are themselves on glycan-covered surface is not.

This searches small rigid-body displacements of the modelled pose and asks for the
smallest move that clears every glycan while keeping as much of the gp120 interface as
possible. The search is rigid and the glycans are rigid, so it answers an accessibility
question, not an affinity one: it tells us whether a glycan-free pose exists nearby, not
what it would be worth.

Scoring per pose:
  glycan clash   any binder heavy atom within 3.0 A of a sugar heavy atom
  self clash     any binder heavy atom within 2.5 A of gp120 protomer A (a pose pushed
                 into the receptor is not a solution)
  retention      binder-gp120 heavy-atom pairs within 4.5 A, against the modelled pose

Run: python glycan_accommodation.py
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

GLYCAN = {"NAG", "BMA", "MAN"}
GP41 = ("B", "F", "J")
VRC01 = ("C", "D", "G", "H", "K", "L")

CLASH_GLY = 3.0        # A, binder vs sugar
CLASH_SELF = 2.5       # A, binder vs its own receptor
CONTACT = 4.5          # A
MAX_T = 4.0            # A, largest translation sampled
MAX_R = 20.0           # degrees, largest rotation sampled
N_POSE = 25000
SEED = 20260211

parser = PDBParser(QUIET=True)
rng = np.random.default_rng(SEED)


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def rot(axis, deg):
    a = axis / np.linalg.norm(axis)
    th = np.deg2rad(deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


tri = next(iter(parser.get_structure("t", os.path.join(SD, "6V8X.pdb"))))
gly, neigh = [], []
for ch in tri:
    for res in ch:
        if res.get_resname().strip() in GLYCAN:
            gly += [a.coord for a in res if a.element != "H"]
        elif res.id[0] == " " and (ch.id in GP41 or (ch.id in ("E", "I"))):
            neigh += [a.coord for a in res if a.element != "H"]
GLY = np.array(gly, dtype=float)
NEIGH = np.array(neigh, dtype=float)
A_xyz = np.array([a.coord for r in tri["A"] if r.id[0] == " "
                  for a in r if a.element != "H"], dtype=float)

t_gly, t_self, t_nb = cKDTree(GLY), cKDTree(A_xyz), cKDTree(NEIGH)

CASES = [("D5", "Latest_Best_split.pdb"), ("D2", "Chai1_Best_split.pdb"),
         ("D2p", "model_4_chai_split.pdb"), ("D3", "Cyclic_Ring_split.pdb")]


def binder(fname):
    st = next(iter(parser.get_structure("d", os.path.join(SD, fname))))
    return np.array([a.coord for r in st["B"] if r.id[0] == " "
                     for a in r if a.element != "H"], dtype=float)


def score(P, self_allowed=CLASH_SELF):
    """(glycan clashes, self clashes, contacts, neighbour clashes) for a pose.

    self_allowed is the receptor approach the designed pose itself achieves. Using a
    fixed 2.5 A instead would hold every decoy to a stricter standard than the pose
    being tested, which biases the search toward poses that sit further off the surface
    and so retain fewer contacts.
    """
    g = int((t_gly.query(P, k=1)[0] < CLASH_GLY).sum())
    s = int((t_self.query(P, k=1)[0] < self_allowed).sum())
    c = int(t_self.query_ball_point(P, CONTACT, return_length=True).sum())
    n = int((t_nb.query(P, k=1)[0] < self_allowed).sum())
    return g, s, c, n


results = {}
for tag, fname in CASES:
    P0 = binder(fname)
    cen = P0.mean(0)
    native_min = float(t_self.query(P0, k=1)[0].min())
    self_allowed = min(CLASH_SELF, native_min)
    g0, s0, c0, n0 = score(P0, self_allowed)
    best = None
    best_minclash = (g0, 0.0, c0)
    # best contacts seen at each glycan-clash count; the single lowest-clash pose is a
    # noisy statistic, while the trade-off curve over a clash budget is stable
    front = {}
    # the modelled pose is itself a candidate: a binder that is already clear needs no move
    if g0 == 0 and s0 == 0 and n0 == 0:
        best = (0.0, -1.0, 0.0, 0.0, c0)

    for _ in range(N_POSE):
        tmag = rng.random() ** 0.5 * MAX_T          # bias toward small moves
        tvec = rng.normal(size=3)
        tvec *= tmag / np.linalg.norm(tvec)
        ang = rng.random() ** 0.5 * MAX_R
        R = rot(rng.normal(size=3), ang)
        P = (P0 - cen) @ R.T + cen + tvec
        g, s, c, n = score(P, self_allowed)
        if g < best_minclash[0] or (g == best_minclash[0] and c > best_minclash[2]):
            best_minclash = (g, tmag, c)
        if s == 0 and n == 0 and c > front.get(g, -1):
            front[g] = c
        if g == 0 and s == 0 and n == 0:
            keep = c / max(c0, 1)
            cand = (tmag + ang / 20.0, -keep, tmag, ang, c)
            if best is None or cand < best:
                best = cand

    # max contacts retained when at most this fraction of the modelled clashes remain
    budgets = {}
    for frac in (0.0, 0.05, 0.10, 0.25, 0.50):
        cap = int(round(frac * g0))
        hits = [cc for gg, cc in front.items() if gg <= cap]
        budgets["%d%%" % round(100 * frac)] = (
            {"max_clashes_allowed": cap,
             "contacts_retained_pct": round(100.0 * max(hits) / max(c0, 1), 1)}
            if hits else {"max_clashes_allowed": cap, "contacts_retained_pct": None})

    row = {"clash_budget_curve": budgets,
           "modelled_pose": {"glycan_clashes": g0, "contacts": c0,
                             "self_clashes": s0, "neighbour_clashes": n0,
                             "closest_receptor_approach_A": round(10 * native_min, 2)},
           "poses_sampled": N_POSE,
           "max_translation_A": MAX_T, "max_rotation_deg": MAX_R}
    if best is None:
        row["glycan_free_pose_found"] = False
        row["lowest_glycan_clashes_reached"] = best_minclash[0]
        row["contacts_at_that_pose"] = best_minclash[2]
        row["contacts_retained_pct"] = round(100.0 * best_minclash[2] / max(c0, 1), 1)
    else:
        _, negkeep, tmag, ang, c = best
        row["glycan_free_pose_found"] = True
        row["translation_A"] = round(tmag, 2)
        row["rotation_deg"] = round(ang, 1)
        row["contacts"] = c
        row["contacts_retained_pct"] = round(-100.0 * negkeep, 1)
    results[tag] = row
    print("   budget curve:", ", ".join(
        "%s->%s" % (k, v["contacts_retained_pct"]) for k, v in budgets.items()))
    print("%-4s modelled: %3d glycan clashes, %4d contacts  ->  %s"
          % (tag, g0, c0,
             ("clear at %.2f A / %.1f deg, %.1f%% contacts kept"
              % (row["translation_A"], row["rotation_deg"], row["contacts_retained_pct"]))
             if row["glycan_free_pose_found"] else
             ("no glycan-free pose within %.0f A; best %d clashes at %.1f%% contacts"
              % (MAX_T, row["lowest_glycan_clashes_reached"],
                 row["contacts_retained_pct"]))))

with io.open(os.path.join(OUT, "glycan_accommodation.json"), "w", encoding="utf-8") as fh:
    json.dump({"search": {"poses": N_POSE, "seed": SEED, "max_translation_A": MAX_T,
                          "max_rotation_deg": MAX_R,
                          "criteria": {"glycan_clash_A": CLASH_GLY,
                                       "self_clash_A": CLASH_SELF,
                                       "contact_A": CONTACT}},
               "note": "rigid binder, rigid glycans, glycans as resolved in 6V8X; "
                       "an accessibility test, not an affinity calculation",
               "partners": results}, fh, indent=2)
print("\nwrote glycan_accommodation.json")
