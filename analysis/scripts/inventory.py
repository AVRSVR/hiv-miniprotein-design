"""Characterize every candidate PDB: chains, residue counts, sequences, ligand/receptor detection."""
import os, glob, json, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

DOWNLOADS = r"C:\Users\anant\Downloads"
PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
OUT = os.path.join(PROJECT, "paper_analysis", "00_inventory")

# gp120 is ~480 aa; designed miniproteins are ~20-120 aa
GP120_MIN_LEN = 300

parser = PDBParser(QUIET=True)


def seq_of(chain):
    s = []
    for res in chain:
        if res.id[0] != " ":
            continue
        try:
            s.append(protein_letters_3to1[res.resname])
        except KeyError:
            s.append("X")
    return "".join(s)


def describe(path):
    try:
        st = parser.get_structure("s", path)
    except Exception as e:
        return {"file": os.path.basename(path), "error": str(e)}
    model = next(iter(st))
    chains = []
    for ch in model:
        sq = seq_of(ch)
        if len(sq) == 0:
            continue
        chains.append({"chain": ch.id, "len": len(sq), "seq": sq})
    if not chains:
        return {"file": os.path.basename(path), "error": "no protein chains"}
    has_receptor = any(c["len"] >= GP120_MIN_LEN for c in chains)
    binders = [c for c in chains if c["len"] < GP120_MIN_LEN]
    return {
        "file": os.path.basename(path),
        "path": path,
        "n_chains": len(chains),
        "chains": chains,
        "type": "complex" if (has_receptor and binders) else ("receptor_only" if has_receptor else "binder_only"),
        "binder_len": binders[0]["len"] if binders else None,
        "binder_seq": binders[0]["seq"] if binders else None,
    }


# Candidate design files worth characterizing (skip unrelated downloads)
NAME_HINTS = ["strategy", "god_particle", "chai", "cyclic", "bispecific", "candidate",
              "mini_cd4", "rosetta_003", "ai_evolved", "ai_true", "complex", "design",
              "latest_best", "s310", "champion", "model_"]

files = []
for p in glob.glob(os.path.join(DOWNLOADS, "*.pdb")):
    n = os.path.basename(p).lower()
    if any(h in n for h in NAME_HINTS):
        files.append(p)
for p in glob.glob(os.path.join(PROJECT, "*.pdb")):
    files.append(p)

results = [describe(p) for p in sorted(set(files))]
results = [r for r in results if "error" not in r]

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "inventory.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"{'FILE':<70} {'TYPE':<14} {'CHAINS':<8} {'BINDER_LEN'}")
print("-" * 110)
for r in sorted(results, key=lambda x: (x["type"], x["file"])):
    print(f"{r['file'][:69]:<70} {r['type']:<14} {r['n_chains']:<8} {r['binder_len']}")
print(f"\nTotal: {len(results)} structures")
