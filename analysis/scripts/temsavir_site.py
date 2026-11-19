"""
Where does the design sit relative to the one approved small molecule that binds gp120?

Fostemsavir is the only licensed HIV-1 attachment inhibitor. Its active form,
temsavir (BMS-626529), binds a pocket under the beta20-beta21 loop and holds Env in the
closed state so CD4 cannot engage it. That pocket is close to the surface D5 occupies,
which raises two questions a reader will have: is the design simply rediscovering the
temsavir site, and could the two be combined?

5U7O is BG505 SOSIP.664 with temsavir bound, the same isolate as the design target, so
its gp120 superposes on 6V8X chain A cleanly and the drug can be carried into the frame
the designs already occupy. We then measure the distance from each design to the drug
and compare the pocket residues with each footprint.

Run: python temsavir_site.py
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

LIG = "83J"              # temsavir / BMS-626529
POCKET = 5.0             # A, gp120 residues lining the drug
FOOTPRINT = 5.0
DESIGNS = [("D5", "Latest_Best_split.pdb"), ("D2", "Chai1_Best_split.pdb"),
           ("D2p", "model_4_chai_split.pdb"), ("D3", "Cyclic_Ring_split.pdb")]

parser = PDBParser(QUIET=True)


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def kabsch(mob, ref):
    mc, rc = mob.mean(0), ref.mean(0)
    V, S, Wt = np.linalg.svd((mob - mc).T @ (ref - rc))
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, rc - (mc @ R)


tri = next(iter(parser.get_structure("t", os.path.join(SD, "6V8X.pdb"))))["A"]
A_res = {r.id[1]: r for r in tri if r.id[0] == " "}
A_ca = {n: r["CA"].coord for n, r in A_res.items() if "CA" in r}

st = next(iter(parser.get_structure("d", os.path.join(SD, "5U7O.pdb"))))
g = st["G"]
drug = [r for r in g if r.get_resname().strip() == LIG]
own_ca = {r.id[1]: r["CA"].coord for r in g if r.id[0] == " " and "CA" in r}
common = sorted(set(own_ca) & set(A_ca))
R, t = kabsch(np.array([own_ca[n] for n in common], dtype=float),
              np.array([A_ca[n] for n in common], dtype=float))
rms = float(np.sqrt(((((np.array([own_ca[n] for n in common]) @ R + t))
                      - np.array([A_ca[n] for n in common])) ** 2).sum(-1).mean()))
print("5U7O gp120 onto 6V8X chain A: %d common CA, fit %.2f A, %d copies of %s"
      % (len(common), rms, len(drug), LIG))

D = np.vstack([heavy(r) @ R + t for r in drug])
# pocket residues, read in the 6V8X frame so they are directly comparable to footprints
pocket = sorted(n for n, r in A_res.items() if heavy(r).size
                and np.linalg.norm(heavy(r)[:, None, :] - D[None, :, :], axis=-1).min()
                < POCKET)
print("pocket residues:", ", ".join("%s%d" % (A_res[n].get_resname().title(), n)
                                    for n in pocket))

rows = {}
for tag, fname in DESIGNS:
    sd = next(iter(parser.get_structure("b", os.path.join(SD, fname))))
    B = np.array([a.coord for r in sd["B"] if r.id[0] == " "
                  for a in r if a.element != "H"], dtype=float)
    fp = sorted(n for n, r in A_res.items() if heavy(r).size
                and np.linalg.norm(heavy(r)[:, None, :] - B[None, :, :], axis=-1).min()
                < FOOTPRINT)
    shared = sorted(set(fp) & set(pocket))
    dmin = float(np.linalg.norm(B[:, None, :] - D[None, :, :], axis=-1).min())
    rows[tag] = {"footprint_residues": len(fp),
                 "min_dist_to_temsavir_A": round(dmin, 2),
                 "clashes_with_temsavir_3A": int((np.linalg.norm(
                     B[:, None, :] - D[None, :, :], axis=-1) < 3.0).sum()),
                 "pocket_residues_shared": ["%s%d" % (A_res[n].get_resname().title(), n)
                                            for n in shared],
                 "pocket_overlap_pct": round(100.0 * len(shared) / len(pocket), 1)}
    print("%-4s closest approach %6.2f A | shares %d of %d pocket residues (%s)"
          % (tag, dmin, len(shared), len(pocket),
             ", ".join(rows[tag]["pocket_residues_shared"]) or "none"))

with io.open(os.path.join(OUT, "temsavir_site.json"), "w", encoding="utf-8") as fh:
    json.dump({"reference": "5U7O, BG505 SOSIP.664 with temsavir (BMS-626529) bound",
               "ligand": LIG, "fit_to_6V8X_A": round(rms, 2), "common_CA": len(common),
               "pocket_residues": ["%s%d" % (A_res[n].get_resname().title(), n)
                                   for n in pocket],
               "designs": rows}, fh, indent=2)
print("\nwrote temsavir_site.json")
