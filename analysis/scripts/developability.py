"""
Developability profile of each designed binder.

A binder that engages its target is not yet a candidate: it also has to express,
stay soluble, and survive storage without chemical degradation. These are the
standard sequence- and structure-level flags used to triage designs before
anything is made.

Sequence level
  MW, isoelectric point, net charge at pH 7.4, GRAVY hydropathy, aromaticity,
  instability index (Guruprasad; >40 predicts a short in vivo half-life)

Chemical liabilities
  unpaired cysteines (disulfide scrambling, covalent aggregation)
  N-glycosylation sequons N-X-[ST], X != P
  NG / NS motifs (asparagine deamidation)
  DG / DP motifs (aspartate isomerisation)
  surface-exposed methionine (oxidation)

Structure level
  secondary-structure content (DSSP), radius of gyration,
  exposed hydrophobic surface and the largest contiguous hydrophobic patch -
  the main structural correlate of aggregation propensity
"""
import os, re, json, warnings
import numpy as np
warnings.filterwarnings("ignore")

import mdtraj as md
from Bio.SeqUtils.ProtParam import ProteinAnalysis

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SCENE = os.path.join(PROJECT, "paper_analysis", "04_figures", "scene")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")

LABEL = {"S310": "D1 (CD4-mimetic control)", "Chai1_Best": "D2 (de novo bispecific)",
         "Cyclic_Ring": "D3 (CD4-mimetic cyclic)", "model_4_chai": "D4 (de novo)",
         "Latest_Best": "D5 (de novo, coreceptor site)"}
HYDROPHOBIC = set("AVLIMFWYC")
EXPOSED_FRAC = 0.25        # relative SASA above which a residue counts as exposed


def liabilities(seq):
    return {
        "n_cys": seq.count("C"),
        "n_free_cys_if_unpaired": seq.count("C") % 2,
        "n_glycosylation_sequons": len(re.findall(r"N[^P][ST]", seq)),
        "n_deamidation_NG_NS": len(re.findall(r"N[GS]", seq)),
        "n_isomerisation_DG_DP": len(re.findall(r"D[GP]", seq)),
        "n_met": seq.count("M"),
    }


def structure_metrics(pdb):
    t = md.load(pdb)
    dssp = md.compute_dssp(t, simplified=True)[0]
    n = len(dssp)
    helix = float((dssp == "H").sum()) / n
    sheet = float((dssp == "E").sum()) / n
    coil = float((dssp == "C").sum()) / n
    rg = float(md.compute_rg(t)[0]) * 10.0

    sasa_res = md.shrake_rupley(t, mode="residue")[0] * 100.0    # A^2
    # Tien et al. maximum accessible surface areas, for relative SASA
    maxacc = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
              "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197,
              "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
              "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}
    res = list(t.topology.residues)
    rel, exposed_phobic, exposed_met = [], 0.0, 0
    patch_run, best_patch = 0, 0
    for i, r in enumerate(res):
        m = maxacc.get(r.name, 200)
        f = sasa_res[i] / m
        rel.append(f)
        code = r.code if r.code else "X"
        if code in HYDROPHOBIC and f > EXPOSED_FRAC:
            exposed_phobic += sasa_res[i]
            patch_run += 1
            best_patch = max(best_patch, patch_run)
        else:
            patch_run = 0
        if r.name == "MET" and f > EXPOSED_FRAC:
            exposed_met += 1

    buried_phobic = sum(1 for i, r in enumerate(res)
                        if (r.code or "X") in HYDROPHOBIC and rel[i] <= EXPOSED_FRAC)
    n_phobic = sum(1 for r in res if (r.code or "X") in HYDROPHOBIC)
    return {
        "helix_frac": round(helix, 3), "sheet_frac": round(sheet, 3),
        "coil_frac": round(coil, 3), "radius_of_gyration_A": round(rg, 2),
        "total_SASA_A2": round(float(sasa_res.sum()), 1),
        "exposed_hydrophobic_SASA_A2": round(float(exposed_phobic), 1),
        "exposed_hydrophobic_frac": round(float(exposed_phobic / max(sasa_res.sum(), 1)), 3),
        "largest_hydrophobic_patch_res": int(best_patch),
        "buried_hydrophobic_frac": round(buried_phobic / max(n_phobic, 1), 3),
        "exposed_methionines": exposed_met,
    }


rows = []
for tag, label in LABEL.items():
    pdb = os.path.join(SCENE, f"binder_{tag}.pdb")
    if not os.path.exists(pdb):
        continue
    t = md.load(pdb)
    seq = "".join((r.code or "X") for r in t.topology.residues)
    pa = ProteinAnalysis(seq.replace("X", "A"))

    row = {"design": label, "tag": tag, "length": len(seq), "sequence": seq,
           "MW_kDa": round(pa.molecular_weight() / 1000.0, 2),
           "pI": round(pa.isoelectric_point(), 2),
           "net_charge_pH7.4": round(pa.charge_at_pH(7.4), 1),
           "GRAVY": round(pa.gravy(), 3),
           "aromaticity": round(pa.aromaticity(), 3),
           "instability_index": round(pa.instability_index(), 1)}
    row["stable_by_instability_index"] = bool(row["instability_index"] < 40)
    row.update(liabilities(seq))
    row.update(structure_metrics(pdb))
    rows.append(row)

with open(os.path.join(OUT, "developability.json"), "w") as f:
    json.dump(rows, f, indent=2)

hdr = (f"{'design':<30}{'len':>5}{'pI':>7}{'charge':>8}{'GRAVY':>8}{'instab':>8}"
       f"{'expo_phob':>11}{'patch':>7}{'glyc':>6}{'Cys':>5}")
print(hdr)
print("-" * len(hdr))
for r in rows:
    print(f"{r['design'][:29]:<30}{r['length']:>5}{r['pI']:>7.2f}"
          f"{r['net_charge_pH7.4']:>8.1f}{r['GRAVY']:>8.2f}{r['instability_index']:>8.1f}"
          f"{r['exposed_hydrophobic_frac']:>11.3f}"
          f"{r['largest_hydrophobic_patch_res']:>7}"
          f"{r['n_glycosylation_sequons']:>6}{r['n_cys']:>5}")
print(f"\nwrote {os.path.join(OUT, 'developability.json')}")
