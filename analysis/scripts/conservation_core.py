"""
Per-position gp120 conservation in HXB2 numbering, computed once and cached.

Both conservation.py and escape_risk.py need this, and until now each carried its
own copy of the alignment code, which is two chances to drift apart. The alignment
lives here instead, and the result is cached to JSON so the pairwise alignments run
once rather than once per design.

Sampling
--------
The original analysis used 47 reviewed Env entries. Reviewed coverage of HIV-1 Env
is tiny, so this uses every UniProt entry for HIV-1 (taxonomy 11676) gene env of
length 700-1000, i.e. full-length gp160 records that actually contain the whole
gp120 region - roughly a hundredfold more sequences, at the cost of curation.

UniProt's env records are not a balanced sample of circulating HIV-1. They are
whatever has been deposited, so closely related sequences from single studies or
transmission chains are over-represented, and pooling them measures conservation in
the database rather than in the epidemic. Subtype balancing is not available as a
fix: 4,439 of the 4,497 records carry no subtype in their organism name, so any
subtype label would be invented rather than read.

What is available is the redundancy itself. Every retained sequence is aligned to the
same 481 HXB2 columns, so identity between any two is a vector comparison. Sequences
are greedy-clustered at 98% identity over those columns and conservation is recomputed
over one representative per cluster, which stops a heavily sequenced lineage from
counting hundreds of times. Both estimates are cached - "pooled" over all retained
sequences and "clustered" over representatives. The clustered figure is the one to
quote.

Group N and O records are dropped: they sit outside group M, the pandemic group and
the only relevant population for an escape argument about circulating virus.

Filtering, in order: group M only, exact duplicates collapsed, under 2% ambiguous
characters, and alignment coverage of the HXB2 gp120 reference of at least 80%.

Run directly to rebuild the cache.
"""
import io
import json
import os
import re
import warnings
from collections import Counter

import numpy as np

warnings.filterwarnings("ignore")
from Bio import SeqIO, Align

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SDIR = os.path.join(PA, "01_structures")
OUT = os.path.join(PA, "05_tables")
CACHE = os.path.join(OUT, "conservation_by_position.json")
MATRIX_CACHE = os.path.join(OUT, "conservation_matrix.npz")
FASTA = os.path.join(SDIR, "hiv_env_full.fasta")

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
POSITIONS = list(range(HXB2_START, HXB2_START + len(HXB2_GP120)))

MIN_COVERAGE = 0.80
MAX_AMBIGUOUS = 0.02
MIN_COL_N = 10          # positions scored in fewer sequences are not reported
CLUSTER_ID = 0.98       # greedy clustering threshold over the aligned HXB2 columns
AMBIG = set("XBZJUO")
GAP = 0                 # sentinel in the aligned residue matrix

_NON_M = re.compile(r"\bgroup\s+[NO]\b")


def is_group_m(description):
    """False for the handful of group N and O records, which are not pandemic HIV-1."""
    return _NON_M.search(description) is None


def _aligner():
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    return al


def _modal(col):
    c = Counter(col)
    return c.most_common(1)[0][1] / sum(c.values())


def _entropy(col):
    c = Counter(col)
    n = sum(c.values())
    return -sum((v / n) * np.log2(v / n) for v in c.values())


def _cluster(M, threshold=CLUSTER_ID):
    """
    Greedy clustering of the aligned residue matrix; returns representative indices.

    Identity between two sequences is counted over the columns where both have a
    residue, so a sequence missing a stretch is not penalised for the gap. The first
    unassigned sequence seeds each cluster, which is the usual greedy scheme.
    """
    n = M.shape[0]
    taken = np.zeros(n, bool)
    reps, sizes = [], []
    for i in range(n):
        if taken[i]:
            continue
        reps.append(i)
        taken[i] = True
        rest = np.where(~taken)[0]
        if rest.size:
            both = (M[i] != GAP) & (M[rest] != GAP)
            same = (M[rest] == M[i]) & both
            denom = both.sum(axis=1)
            ident = np.where(denom > 0, same.sum(axis=1) / np.maximum(denom, 1), 0.0)
            hit = rest[ident >= threshold]
            taken[hit] = True
            sizes.append(1 + hit.size)
        else:
            sizes.append(1)
    return np.array(reps), np.array(sizes)


def compute(fasta=FASTA, verbose=True):
    """Align every retained sequence to HXB2 gp120 and tabulate each position."""
    recs = list(SeqIO.parse(fasta, "fasta"))
    seen, kept = set(), []
    n_dup = n_ambig = n_nonm = 0
    for r in recs:
        if not is_group_m(r.description):
            n_nonm += 1
            continue
        s = str(r.seq).upper()
        if s in seen:
            n_dup += 1
            continue
        seen.add(s)
        if sum(ch in AMBIG for ch in s) / max(len(s), 1) > MAX_AMBIGUOUS:
            n_ambig += 1
            continue
        kept.append(s)

    al = _aligner()
    idx = {p: j for j, p in enumerate(POSITIONS)}
    rows, n_lowcov = [], 0
    for s in kept:
        try:
            aln = al.align(HXB2_GP120, s)[0]
        except Exception:
            continue
        row = np.zeros(len(POSITIONS), dtype=np.uint8)
        pos, cov = HXB2_START, 0
        for cr, cq in zip(aln[0], aln[1]):
            if cr == "-":
                continue
            if cq != "-" and cq not in AMBIG:
                row[idx[pos]] = ord(cq)
                cov += 1
            pos += 1
        if cov / len(POSITIONS) < MIN_COVERAGE:
            n_lowcov += 1
            continue
        rows.append(row)

    M = np.array(rows, dtype=np.uint8)
    reps, sizes = _cluster(M)
    Mc = M[reps]
    # the clustered matrix itself, for analyses that need per-sequence residues
    # rather than per-position frequencies
    np.savez_compressed(MATRIX_CACHE, Mc=Mc, positions=np.array(POSITIONS),
                        cluster_sizes=sizes)

    per_pos = {}
    for p, j in idx.items():
        col = [chr(c) for c in M[:, j] if c != GAP]
        if len(col) < MIN_COL_N:
            continue
        colc = [chr(c) for c in Mc[:, j] if c != GAP]
        per_pos[p] = {
            "n": len(col),
            "n_clustered": len(colc),
            "pooled": round(_modal(col), 4),
            "clustered": round(_modal(colc), 4) if len(colc) >= MIN_COL_N else None,
            "entropy": round(_entropy(col), 4),
        }

    meta = {
        "fasta": os.path.basename(fasta),
        "records_in_file": len(recs),
        "dropped_group_N_O": n_nonm,
        "duplicates_collapsed": n_dup,
        "dropped_ambiguous": n_ambig,
        "dropped_low_coverage": n_lowcov,
        "sequences_used": int(M.shape[0]),
        "clusters_at_98pct_identity": int(reps.size),
        "largest_cluster": int(sizes.max()) if reps.size else 0,
        "min_coverage": MIN_COVERAGE,
        "cluster_identity": CLUSTER_ID,
        "positions_scored": len(per_pos),
    }
    if verbose:
        print("records %d | group N/O %d | duplicates %d | ambiguous %d | "
              "low coverage %d | used %d"
              % (len(recs), n_nonm, n_dup, n_ambig, n_lowcov, M.shape[0]))
        print("clusters at %.0f%% identity: %d (largest %d sequences)"
              % (100 * CLUSTER_ID, reps.size, sizes.max() if reps.size else 0))
        print("positions scored:", len(per_pos))
    return meta, per_pos


def matrix(rebuild=False):
    """The clustered aligned residue matrix: one row per cluster representative,
    one column per HXB2 gp120 position, ASCII codes with 0 for a gap.

    Per-position conservation collapses the sequence set to a frequency, which
    cannot answer how often a whole set of positions is intact in the same virus.
    That question needs the rows, so they are cached alongside the summary.
    """
    if rebuild or not os.path.exists(MATRIX_CACHE):
        compute()
    d = np.load(MATRIX_CACHE)
    return d["Mc"], {int(p): j for j, p in enumerate(d["positions"])}, d["cluster_sizes"]


def load(prefer="clustered", rebuild=False):
    """{position: conservation} using the cache, building it when absent."""
    if rebuild or not os.path.exists(CACHE):
        meta, per_pos = compute()
        with io.open(CACHE, "w", encoding="utf-8") as fh:
            json.dump({"meta": meta,
                       "per_position": {str(k): v for k, v in per_pos.items()}},
                      fh, indent=1)
    with io.open(CACHE, encoding="utf-8") as fh:
        d = json.load(fh)
    out = {}
    for k, v in d["per_position"].items():
        c = v.get(prefer)
        out[int(k)] = v["pooled"] if c is None else c
    return out, d["meta"]


if __name__ == "__main__":
    cons, meta = load(rebuild=True)
    pooled, _ = load(prefer="pooled")
    print("\nmean conservation: pooled %.1f%% | clustered %.1f%% over %d positions"
          % (100 * np.mean(list(pooled.values())),
             100 * np.mean(list(cons.values())), len(cons)))
    print("cache:", CACHE)
