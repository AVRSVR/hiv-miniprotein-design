"""
Are the targeted epitopes actually conserved?

The manuscript's escape argument rests on the epitopes being functionally
constrained. That is testable: align reviewed HIV-1 Env sequences spanning group M
subtypes, map each column onto HXB2 numbering via the HXB2 reference sequence, and
compare per-residue conservation at each design's gp120 contact residues against
gp120 as a whole.

Conservation is reported as the frequency of the modal residue (Shannon entropy is
also given), computed only over sequences with a non-gap character at that column.
"""
import os, json, warnings
from collections import Counter
import numpy as np
warnings.filterwarnings("ignore")
from Bio import SeqIO, Align

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SDIR = os.path.join(PROJECT, "paper_analysis", "01_structures")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")

# Superseded: the 47-sequence reviewed set below is no longer the source of
# conservation. conservation_core reads hiv_env_full.fasta (4,497 full-length
# UniProt records) and caches the per-position result. Kept only so the HXB2
# reference sequence that follows stays next to its provenance.
_UNUSED_FASTA = os.path.join(SDIR, "hiv_env.fasta")

# HXB2 Env (P04578) gp120 region, residues 31-511 in HXB2 numbering.
HXB2_GP120 = (
    "TEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQ"
    "MHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFF"
    "YKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTH"
    "GIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGK"
    "IGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTW"
    "FNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIF"
    "RPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"
)
HXB2_START = 31

LABEL = {"S310": "D1 (CD4-mimetic control)", "Chai1_Best": "D2 (de novo bispecific)",
         "Cyclic_Ring": "D3 (CD4-mimetic cyclic)", "model_4_chai": "D4 (de novo)",
         "Latest_Best": "D5 (de novo, coreceptor site)"}
CD4BS = [368, 370, 371, 425, 427, 429, 457, 473, 280]

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conservation_core

# The alignment itself lives in conservation_core so that this script and
# escape_risk.py cannot drift apart; it is cached after the first build.
CONS, META = conservation_core.load()            # 98%-identity-clustered
POOLED, _ = conservation_core.load(prefer="pooled")
ENTROPY = conservation_core.load(prefer="entropy")[0]
print("%d Env sequences -> %d clusters at %.0f%% identity, %d positions scored"
      % (META["sequences_used"], META["clusters_at_98pct_identity"],
         100 * META["cluster_identity"], META["positions_scored"]))
print("group N/O dropped %d | duplicates %d | ambiguous %d | low coverage %d"
      % (META["dropped_group_N_O"], META["duplicates_collapsed"],
         META["dropped_ambiguous"], META["dropped_low_coverage"]))
per_pos = CONS


def stats(p):
    if p not in CONS:
        return None
    return {"n": META["sequences_used"], "conservation": CONS[p],
            "pooled": POOLED.get(p), "entropy": ENTROPY.get(p)}


allp = [stats(p) for p in per_pos]
allp = [x for x in allp if x]
gp120_mean = float(np.mean([x["conservation"] for x in allp]))
print(f"gp120 overall: mean conservation {100*gp120_mean:.1f}% "
      f"over {len(allp)} positions")

cd4bs_vals = [stats(p)["conservation"] for p in CD4BS if stats(p)]
print(f"CD4bs hotspots: mean conservation {100*np.mean(cd4bs_vals):.1f}% "
      f"({len(cd4bs_vals)} positions)\n")

iface = json.load(open(os.path.join(SDIR, "interface_metrics.json")))
rows = []
print(f"{'design':<32}{'iface res':>10}{'mean cons':>11}{'>=90% cons':>12}{'vs gp120':>10}")
print("-" * 76)
for d in iface:
    if d["design"] == "docked_complex":
        continue
    vals, hi = [], 0
    for p in d["iface_res_receptor"]:
        s = stats(p)
        if s:
            vals.append(s["conservation"])
            hi += s["conservation"] >= 0.90
    if not vals:
        continue
    mean_c = float(np.mean(vals))
    rows.append({"design": LABEL.get(d["design"], d["design"]),
                 "n_iface_res_scored": len(vals),
                 "mean_conservation_pct": round(100 * mean_c, 1),
                 "n_highly_conserved_ge90pct": hi,
                 "pct_highly_conserved": round(100 * hi / len(vals), 1),
                 "delta_vs_gp120_mean_pct": round(100 * (mean_c - gp120_mean), 1)})
    print(f"{LABEL.get(d['design'], d['design'])[:31]:<32}{len(vals):>10}"
          f"{100*mean_c:>10.1f}%{hi:>7}/{len(vals):<4}"
          f"{100*(mean_c-gp120_mean):>+9.1f}")

with open(os.path.join(OUT, "epitope_conservation.json"), "w") as f:
    json.dump({"n_sequences": META["sequences_used"],
               "n_clusters": META["clusters_at_98pct_identity"],
               "estimator": "modal residue frequency, 98%-identity-clustered",
               "gp120_mean_conservation_pct": round(100 * gp120_mean, 1),
               "cd4bs_hotspot_mean_conservation_pct": round(100 * float(np.mean(cd4bs_vals)), 1),
               "per_design": rows}, f, indent=2)
