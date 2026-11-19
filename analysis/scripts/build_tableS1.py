"""
Table S1: cross-strain conservation of each footprint, against the CD4 and 17b epitopes.

The earlier version of this table predated the expansion of the sequence set and scored
only the designs, which left the numbers without a scale. Here every footprint is
defined the same way, as the gp120 residues within 5 A of the partner, and conservation
comes from the shared cache in conservation_core, so the designs and the two
experimentally validated partners are measured identically.

Run: python build_tableS1.py
"""
import io
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_core
from Bio.PDB import PDBParser

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = os.path.join(PA, "01_structures")
OUT = os.path.join(PA, "05_tables")
FOOTPRINT = 5.0

CASES = [("CD4 (experimental reference)", "1GC1.pdb", ["C"], "G"),
         ("17b Fab (experimental reference)", "1GC1.pdb", ["H", "L"], "G"),
         ("D1 (CD4-mimetic control)", "S310_split.pdb", ["B"], "A"),
         ("D2 (de novo bispecific)", "Chai1_Best_split.pdb", ["B"], "A"),
         ("D2' (alternative pose)", "model_4_chai_split.pdb", ["B"], "A"),
         ("D3 (CD4-mimetic cyclic)", "Cyclic_Ring_split.pdb", ["B"], "A"),
         ("D5 (de novo, coreceptor site)", "Latest_Best_split.pdb", ["B"], "A")]

parser = PDBParser(QUIET=True)
cons, meta = conservation_core.load()
print("conservation: %d sequences, %d clusters at %.0f%% identity"
      % (meta["sequences_used"], meta["clusters_at_98pct_identity"],
         100 * meta["cluster_identity"]))


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def footprint(path, binder_chains, rec_chain):
    st = next(iter(parser.get_structure("x", os.path.join(SD, path))))
    B = np.vstack([heavy(r) for c in binder_chains for r in st[c] if r.id[0] == " "])
    g = st[rec_chain]
    return [r.id[1] for r in g if r.id[0] == " " and heavy(r).size
            and np.linalg.norm(heavy(r)[:, None, :] - B[None, :, :], axis=-1).min()
            < FOOTPRINT]


gp120_mean = 100 * float(np.mean(list(cons.values())))
rows = []
for label, path, chains, rec in CASES:
    fp = footprint(path, chains, rec)
    v = [cons[n] for n in fp if n in cons]
    if not v:
        continue
    rows.append({"Partner": label,
                 "Footprint residues": len(fp),
                 "Scored against the alignment": len(v),
                 "Mean conservation (%)": round(100 * float(np.mean(v)), 1),
                 "Residues >= 90% conserved": int(sum(1 for x in v if x >= 0.90)),
                 "Fraction >= 90% conserved (%)": round(
                     100.0 * sum(1 for x in v if x >= 0.90) / len(v), 1),
                 "Difference from gp120 mean (pp)": round(
                     100 * float(np.mean(v)) - gp120_mean, 1)})

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "tableS1_epitope_conservation.csv"), index=False)
with io.open(os.path.join(OUT, "tableS1_epitope_conservation.md"), "w",
             encoding="utf-8") as fh:
    fh.write(df.to_markdown(index=False))
    fh.write("\n\nConservation is the modal-residue frequency over %d cluster "
             "representatives drawn from %d full-length group M Env records; the gp120 "
             "mean over all %d scored positions is %.1f%%. Footprints are the gp120 "
             "residues within %.1f A of the partner, defined identically for the designs "
             "and for the two experimental references.\n"
             % (meta["clusters_at_98pct_identity"], meta["sequences_used"],
                meta["positions_scored"], gp120_mean, FOOTPRINT))
print(df.to_string(index=False))
print("\ngp120 mean %.1f%%" % gp120_mean)
