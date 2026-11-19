"""
Would this binder stick to anything? A counter-screen against unrelated surfaces.

Limitation 3 of this work has been that specificity is untested. 2.13 removed the
mechanism by which an idealised Glu/Lys-rich helical binder would adhere indiscriminately
to charged surfaces, by showing its interface is van der Waals-driven with an
electrostatic share inside the experimental range, but it never asked the question
directly: put the same binder on a protein it was not designed against and see what it
scores.

The test places the binder face-on at randomly chosen surface points of each target,
minimises every placement identically, and scores it with the same MM/GBSA used
throughout. The comparison that matters is between targets, not between a placement and
the designed pose, so gp120 is screened by exactly the same procedure as the decoy
targets and its designed pose is reported separately rather than mixed into its own
distribution.

Targets:
  gp120     the design target, screened like the others
  1AO6      human serum albumin, the most abundant plasma protein and the off-target
            that matters most for an injected agent
  1DPX      hen egg lysozyme, small and strongly basic, the adversarial surface for a
            binder suspected of charge-driven adhesion
  1UBQ      ubiquitin, small and near-neutral, a size and shape control

Limits, stated plainly: N face-on placements per target is not a docking search, so each
target's best score is an upper bound on what exhaustive docking would find, and the
bound is loose for a large target. The procedure is identical across targets, so the
comparison between them is fair even though no single number is an affinity. A binder
that scored as well on albumin as on gp120 would be damning; one that does not has
passed a necessary test, not a sufficient one.

Because gp120 is screened by the same procedure, this also answers a second question
the work had left open: whether the designed pose is better than an arbitrary placement
on its own target, which is the discrimination a real binder must show.

Run: python counter_screen.py [design_tag ...]
"""
import io
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from scipy.spatial import cKDTree
from openmm.unit import nanometer, kilojoule_per_mole

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmgbsa import build_topology, make_context, subset, energy, SPLIT, OUT
CONTACT = 0.45          # nm
MIN_STEPS = 1500        # minimisation steps applied identically to every placement

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRATCH = os.path.join(PA, "06_scratch")
os.makedirs(SCRATCH, exist_ok=True)

KCAL = 1.0 / 4.184
N_PLACE = 8
SEED = 20261106

TARGETS = [("gp120", None), ("HSA_1AO6", "1AO6.pdb"), ("lysozyme_1DPX", "1DPX.pdb"),
           ("ubiquitin_1UBQ", "1UBQ.pdb")]
DESIGNS = {"Latest_Best": "D5 de novo, coreceptor site"}

rng = np.random.default_rng(SEED)


def align(a, b):
    """Rotation taking unit vector a onto unit vector b."""
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * (1.0 / (1.0 + c))


def spin(axis, deg):
    a = axis / np.linalg.norm(axis)
    th = np.deg2rad(deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def place_face_on(L0, face_mask, face_dir, R, tree_r, min_allowed,
                  min_contacts, rng):
    """One decoy: present the binder's own binding face to a random receptor surface
    point, spin it freely about that approach, and slide in until first contact.

    Approaching along a random direction with a random orientation almost always makes
    first contact on a single protruding atom, which produces decoys with too few
    contacts to be comparable. Presenting the binding face keeps the decoys face-on, so
    they are placements a binder could plausibly adopt rather than glancing touches.
    """
    cen_r = R.mean(0)
    for _ in range(40):
        t = R[rng.integers(len(R))]
        n = t - cen_r
        nrm = np.linalg.norm(n)
        if nrm < 1e-6:
            continue
        n /= nrm
        M = spin(n, rng.uniform(0, 360)) @ align(face_dir, -n)
        body = (L0 - L0.mean(0)) @ M.T
        face = ((L0[face_mask] - L0.mean(0)) @ M.T).mean(0)
        base = t + n * 0.8 - face
        for step in np.arange(0.0, 2.0, 0.04):
            P = body + base - n * step
            near = tree_r.query(P, k=1)[0]
            if near.min() < min_allowed:
                P = body + base - n * max(step - 0.04, 0.0)
                near = tree_r.query(P, k=1)[0]
                nc = int((near < CONTACT).sum())
                if near.min() >= min_allowed and nc >= min_contacts:
                    return P, nc
                break
    return None, 0


def first_protein_chain(path):
    """Heavy-atom records of the first protein chain, waters and heteroatoms dropped."""
    keep, chain = [], None
    for line in io.open(path, encoding="utf-8", errors="ignore"):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM"):
            continue
        c = line[21]
        if chain is None:
            chain = c
        if c != chain:
            continue
        if line[76:78].strip() == "H":
            continue
        keep.append(line)
    return keep


def merged_pdb(target_lines, binder_lines, out_path):
    with io.open(out_path, "w", encoding="utf-8") as fh:
        for ln in target_lines:
            fh.write(ln[:21] + "A" + ln[22:])
        fh.write("TER\n")
        for ln in binder_lines:
            fh.write(ln[:21] + "B" + ln[22:])
        fh.write("TER\nEND\n")


def binder_lines_of(tag):
    out, started = [], False
    for line in io.open(os.path.join(SPLIT, f"{tag}_split.pdb"), encoding="utf-8"):
        if not line.startswith("ATOM"):
            continue
        if line[21] == "B":
            started = True
            if line[76:78].strip() != "H":
                out.append(line)
        elif started:
            break
    return out


def save(designs, partial=False):
    with io.open(os.path.join(OUT, "counter_screen.json"), "w", encoding="utf-8") as fh:
        json.dump({"placements_per_target": N_PLACE, "seed": SEED,
                   "minimisation_steps": MIN_STEPS, "complete": not partial,
                   "method": "binder presented face-on at random target surface points, "
                             "slid to contact, minimised identically, scored by MM/GBSA; "
                             "identical procedure on every target including gp120",
                   "limits": "N placements is not a docking search, so each best score "
                             "is an upper bound; the comparison between targets is fair "
                             "because the procedure is identical",
                   "designs": designs}, fh, indent=2)


def screen(tag, label):
    rows = {}
    designed = None
    for tname, tfile in TARGETS:
        if tfile is None:
            merged = os.path.join(SPLIT, f"{tag}_split.pdb")
        else:
            merged = os.path.join(SCRATCH, f"cs_{tag}_{tname}.pdb")
            merged_pdb(first_protein_chain(os.path.join(SCRATCH, tfile)),
                       binder_lines_of(tag), merged)
        top, pos = build_topology(merged)
        atoms = list(top.atoms())
        rec = [i for i, a in enumerate(atoms) if a.residue.chain.id == "A"]
        lig = [i for i, a in enumerate(atoms) if a.residue.chain.id == "B"]
        hr = [i for i in rec if atoms[i].element.symbol != "H"]
        hl = [i for i in lig if atoms[i].element.symbol != "H"]
        X = np.array(pos.value_in_unit(nanometer))

        top_r, pos_r = subset(top, pos, "A")
        top_l, pos_l = subset(top, pos, "B")
        sim_c = make_context(top, pos)
        sim_r = make_context(top_r, pos_r)
        sim_l = make_context(top_l, pos_l)

        def relax_and_score(full_xyz):
            sim_c.context.setPositions(full_xyz)
            sim_c.minimizeEnergy(maxIterations=MIN_STEPS)
            st = sim_c.context.getState(getPositions=True, getEnergy=True)
            xyz = np.array(st.getPositions().value_in_unit(nanometer))
            e_c = st.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            return (e_c - energy(sim_r, xyz[rec]) - energy(sim_l, xyz[lig])) * KCAL, xyz

        if tfile is None:
            designed, Xmin = relax_and_score(X)
            X = Xmin
            print("  %-16s designed pose %8.1f kcal/mol" % (tname, designed))
        else:
            # start the binder clear of the target so the first placement is not a clash
            X = X.copy()

        R, L0 = X[hr], X[hl]
        tree = cKDTree(R)
        near = tree.query(L0, k=1)[0]
        if tfile is None:
            face = near < CONTACT
            min_allowed = float(near.min())
            min_contacts = max(int(face.sum()) // 3, 15)
            screen.face = face
            screen.min_allowed = min_allowed
            screen.min_contacts = min_contacts
        face, min_allowed = screen.face, screen.min_allowed
        min_contacts = screen.min_contacts
        face_dir = L0[face].mean(0) - L0.mean(0)
        face_dir /= np.linalg.norm(face_dir)

        vals, tries = [], 0
        while len(vals) < N_PLACE and tries < N_PLACE * 14:
            tries += 1
            P, nc = place_face_on(L0, face, face_dir, R, tree, min_allowed,
                                  min_contacts, rng)
            if P is None:
                continue
            A = L0 - L0.mean(0)
            B = P - P.mean(0)
            U, S, Vt = np.linalg.svd(A.T @ B)
            d = np.sign(np.linalg.det(U @ Vt))
            M = U @ np.diag([1.0, 1.0, d]) @ Vt
            full = X.copy()
            full[lig] = X[lig] @ M + (P.mean(0) - L0.mean(0) @ M)
            v, _ = relax_and_score(full)
            if np.isfinite(v) and abs(v) < 1e4:
                vals.append(float(v))
        a = np.array(vals)
        rows[tname] = {"n_placements": len(a), "attempts": tries,
                       "best_kcal": round(float(a.min()), 1) if len(a) else None,
                       "median_kcal": round(float(np.median(a)), 1) if len(a) else None,
                       "mean_kcal": round(float(a.mean()), 1) if len(a) else None,
                       "sd_kcal": round(float(a.std(ddof=1)), 1) if len(a) > 1 else None}
        print("  %-16s %2d placements: best %8.1f  median %8.1f"
              % (tname, len(a), a.min() if len(a) else float("nan"),
                 np.median(a) if len(a) else float("nan")))
        save({tag: {"label": label,
                    "designed_pose_on_gp120_kcal": round(float(designed), 1),
                    "targets": rows}}, partial=True)
    return {"label": label, "designed_pose_on_gp120_kcal": round(float(designed), 1),
            "targets": rows}


if __name__ == "__main__":
    tags = sys.argv[1:] or list(DESIGNS)
    out = {}
    for t in tags:
        print(t)
        out[t] = screen(t, DESIGNS.get(t, t))
    save(out)
    print("\nwrote counter_screen.json")
