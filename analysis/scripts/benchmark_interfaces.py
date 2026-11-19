"""
Calibrate the designed interfaces against experimental ones.

Interface metrics only mean something relative to a reference. The same code that
characterized the design complexes is run here on two experimentally determined
gp120 interfaces from PDB 1GC1 - gp120/CD4 and gp120/17b Fab - plus the gp120/b12
complex (2NY7), an antibody that binds the CD4-binding site. Designed interfaces
that fall inside the range spanned by these are physically plausible; ones far
outside are not.
"""
import os, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SDIR = os.path.join(PROJECT, "paper_analysis", "01_structures")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")

CONTACT_CUT, HBOND_CUT, SALT_CUT = 4.5, 3.5, 4.0
POS = {"ARG": ["NE", "NH1", "NH2"], "LYS": ["NZ"], "HIS": ["ND1", "NE2"]}
NEG = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}

parser = PDBParser(QUIET=True)
sr = ShrakeRupley()


def heavy(chains, model):
    out = []
    for c in chains:
        if c in model:
            out += [a for a in model[c].get_atoms()
                    if a.element != "H" and a.get_parent().id[0] == " "]
    return out


def sasa_of(path, chains):
    st = parser.get_structure("x", path)
    m = next(iter(st))
    for ch in list(m):
        if ch.id not in chains:
            m.detach_child(ch.id)
    for ch in list(m):
        for r in list(ch):
            if r.id[0] != " ":
                ch.detach_child(r.id)
    sr.compute(m, level="C")
    return sum(c.sasa for c in m)


def analyze(path, rec_chains, lig_chains, name):
    st = parser.get_structure(name, path)
    m = next(iter(st))
    R, L = heavy(rec_chains, m), heavy(lig_chains, m)
    if not R or not L:
        return None
    rc = np.array([a.coord for a in R])
    lc = np.array([a.coord for a in L])
    d = np.linalg.norm(rc[:, None, :] - lc[None, :, :], axis=-1)
    pairs = np.argwhere(d < CONTACT_CUT)

    nhb = sum(1 for i, j in pairs if d[i, j] < HBOND_CUT
              and R[i].element in ("N", "O") and L[j].element in ("N", "O"))
    nsalt = 0
    for i, j in pairs:
        if d[i, j] > SALT_CUT:
            continue
        rr, ll = R[i].get_parent().resname, L[j].get_parent().resname
        rp = R[i].name in POS.get(rr, []); rn = R[i].name in NEG.get(rr, [])
        lp = L[j].name in POS.get(ll, []); ln = L[j].name in NEG.get(ll, [])
        if (rp and ln) or (rn and lp):
            nsalt += 1

    dsasa = (sasa_of(path, rec_chains) + sasa_of(path, lig_chains)
             - sasa_of(path, list(rec_chains) + list(lig_chains)))

    nres_r = len({(R[i].get_parent().get_parent().id, R[i].get_parent().id[1])
                  for i, _ in pairs})
    nres_l = len({(L[j].get_parent().get_parent().id, L[j].get_parent().id[1])
                  for _, j in pairs})
    return {"interface": name, "n_contacts": int(len(pairs)), "n_hbonds": int(nhb),
            "n_saltbridges": int(nsalt), "dSASA_A2": round(float(dsasa), 1),
            "n_iface_res_gp120": nres_r, "n_iface_res_partner": nres_l}


refs = []
g1 = os.path.join(SDIR, "1GC1.pdb")
if os.path.exists(g1):
    for chains, nm in ((("C",), "gp120 / CD4 D1D2 (1GC1, experimental)"),
                       (("L", "H"), "gp120 / 17b Fab (1GC1, experimental)")):
        r = analyze(g1, ("G",), chains, nm)
        if r:
            refs.append(r)

# b12: CD4bs-directed broadly neutralizing antibody
b12 = os.path.join(SDIR, "2NY7.pdb")
if not os.path.exists(b12):
    os.system(f'curl -s -o "{b12}" https://files.rcsb.org/download/2NY7.pdb')
if os.path.exists(b12) and os.path.getsize(b12) > 10000:
    r = analyze(b12, ("G",), ("L", "H"), "gp120 / b12 Fab (2NY7, experimental)")
    if r:
        refs.append(r)

designs = json.load(open(os.path.join(SDIR, "interface_metrics.json")))
LABEL = {"S310": "D1 (CD4-mimetic control)", "Chai1_Best": "D2 (de novo bispecific)",
         "Cyclic_Ring": "D3 (CD4-mimetic cyclic)", "model_4_chai": "D4 (de novo)",
         "Latest_Best": "D5 (de novo, coreceptor site)"}

rows = []
for r in refs:
    rows.append({**r, "type": "experimental"})
for d in designs:
    if d["design"] == "docked_complex":
        continue
    rows.append({"interface": LABEL.get(d["design"], d["design"]), "type": "designed",
                 "n_contacts": d["n_contacts"], "n_hbonds": d["n_hbonds"],
                 "n_saltbridges": d["n_saltbridges"], "dSASA_A2": d["dSASA_A2"],
                 "n_iface_res_gp120": d["n_iface_res_receptor"],
                 "n_iface_res_partner": d["n_iface_res_binder"]})

with open(os.path.join(OUT, "interface_benchmark.json"), "w") as f:
    json.dump(rows, f, indent=2)

exp = [r for r in rows if r["type"] == "experimental"]
if exp:
    lo = min(r["dSASA_A2"] for r in exp); hi = max(r["dSASA_A2"] for r in exp)
    print(f"experimental BSA range: {lo:.0f} - {hi:.0f} A^2\n")

print(f"{'interface':<46}{'type':<14}{'contacts':>9}{'H-bond':>8}{'salt':>6}{'BSA':>9}")
print("-" * 94)
for r in sorted(rows, key=lambda x: (x["type"] != "experimental", -x["dSASA_A2"])):
    print(f"{r['interface'][:45]:<46}{r['type']:<14}{r['n_contacts']:>9}"
          f"{r['n_hbonds']:>8}{r['n_saltbridges']:>6}{r['dSASA_A2']:>9.0f}")

if exp:
    print()
    for r in rows:
        if r["type"] == "designed":
            inside = lo * 0.6 <= r["dSASA_A2"] <= hi * 1.4
            print(f"  {r['interface'][:45]:<46} BSA {'within' if inside else 'OUTSIDE'} "
                  f"experimental range")
