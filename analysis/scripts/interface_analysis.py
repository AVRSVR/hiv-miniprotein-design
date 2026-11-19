"""
Interface characterization of designed miniprotein / gp120 complexes.

Complex PDBs carry receptor (gp120, HXB2 numbering) and the designed binder
in the same chain ID, separated by TER records. This splits them at the last
TER, then quantifies the interface:

  - buried solvent-accessible surface area (dSASA)
  - heavy-atom contacts (< 4.5 A)
  - hydrogen bonds (N/O donor-acceptor < 3.5 A)
  - salt bridges (charged groups < 4.0 A)
  - CD4-binding-site hotspot engagement (which known hotspots the design touches)
"""
import os, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
DOWNLOADS = r"C:\Users\anant\Downloads"
OUT = os.path.join(PROJECT, "paper_analysis", "01_structures")
os.makedirs(OUT, exist_ok=True)

CONTACT_CUT = 4.5
HBOND_CUT = 3.5
SALT_CUT = 4.0

# CD4-binding-site hotspots (HXB2 numbering)
HOTSPOTS = [368, 370, 427, 457, 429, 421, 473, 461, 257, 280, 255, 371, 105, 425]
HOTSPOT_LABELS = {368: "D368", 370: "E370", 427: "W427", 457: "D457", 429: "R429",
                  421: "K421", 473: "G473", 461: "E461", 257: "T257", 280: "N280",
                  255: "V255", 371: "I371", 105: "Q105", 425: "N425"}

POS_RES = {"ARG": ["NE", "NH1", "NH2"], "LYS": ["NZ"], "HIS": ["ND1", "NE2"]}
NEG_RES = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}

parser = PDBParser(QUIET=True)
sr = ShrakeRupley()


def split_complex(path, out_dir, tag):
    """Split a receptor+ligand PDB at the last TER into two chain-labeled files."""
    lines = open(path, errors="ignore").readlines()
    ter_idx = [i for i, l in enumerate(lines) if l.startswith("TER")]
    atom_idx = [i for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))]
    if not ter_idx or not atom_idx:
        return None
    # last TER that still has atoms after it
    split_at = None
    for t in reversed(ter_idx):
        if any(a > t for a in atom_idx):
            split_at = t
            break
    if split_at is None:
        return None

    rec, lig = [], []
    for i, l in enumerate(lines):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[17:20].strip() in ("HOH", "WAT", "NA", "CL"):
            continue
        if i < split_at:
            rec.append(l[:21] + "A" + l[22:])
        else:
            lig.append(l[:21] + "B" + l[22:])
    if not rec or not lig:
        return None

    combo = os.path.join(out_dir, f"{tag}_split.pdb")
    with open(combo, "w") as f:
        f.writelines(rec)
        f.write("TER\n")
        f.writelines(lig)
        f.write("TER\nEND\n")
    return combo


def residue_key(res):
    return (res.get_parent().id, res.id[1], res.resname)


def analyze(path, tag):
    st = parser.get_structure(tag, path)
    model = next(iter(st))
    if "A" not in model or "B" not in model:
        return None
    rec, lig = model["A"], model["B"]

    rec_atoms = [a for a in rec.get_atoms() if a.element != "H"]
    lig_atoms = [a for a in lig.get_atoms() if a.element != "H"]
    if not rec_atoms or not lig_atoms:
        return None

    rc = np.array([a.coord for a in rec_atoms])
    lc = np.array([a.coord for a in lig_atoms])
    d = np.linalg.norm(rc[:, None, :] - lc[None, :, :], axis=-1)

    # contacts
    pairs = np.argwhere(d < CONTACT_CUT)
    n_contacts = len(pairs)

    rec_iface, lig_iface = set(), set()
    for i, j in pairs:
        rec_iface.add(residue_key(rec_atoms[i].get_parent()))
        lig_iface.add(residue_key(lig_atoms[j].get_parent()))

    # hydrogen bonds (N/O - N/O within cutoff)
    n_hb = 0
    for i, j in pairs:
        if d[i, j] < HBOND_CUT:
            if rec_atoms[i].element in ("N", "O") and lig_atoms[j].element in ("N", "O"):
                n_hb += 1

    # salt bridges
    n_salt = 0
    for i, j in pairs:
        if d[i, j] > SALT_CUT:
            continue
        ra, la = rec_atoms[i], lig_atoms[j]
        rr, ll = ra.get_parent().resname, la.get_parent().resname
        r_pos = rr in POS_RES and ra.name in POS_RES.get(rr, [])
        r_neg = rr in NEG_RES and ra.name in NEG_RES.get(rr, [])
        l_pos = ll in POS_RES and la.name in POS_RES.get(ll, [])
        l_neg = ll in NEG_RES and la.name in NEG_RES.get(ll, [])
        if (r_pos and l_neg) or (r_neg and l_pos):
            n_salt += 1

    # buried SASA
    st_c = parser.get_structure("c", path)
    mc = next(iter(st_c))
    sr.compute(mc, level="C")
    sasa_complex = sum(c.sasa for c in mc)

    st_r = parser.get_structure("r", path)
    mr = next(iter(st_r))
    for ch in list(mr):
        if ch.id != "A":
            mr.detach_child(ch.id)
    sr.compute(mr, level="C")
    sasa_rec = sum(c.sasa for c in mr)

    st_l = parser.get_structure("l", path)
    ml = next(iter(st_l))
    for ch in list(ml):
        if ch.id != "B":
            ml.detach_child(ch.id)
    sr.compute(ml, level="C")
    sasa_lig = sum(c.sasa for c in ml)

    dsasa = sasa_rec + sasa_lig - sasa_complex

    # hotspot engagement (receptor side, HXB2 numbering)
    engaged = sorted({r[1] for r in rec_iface if r[1] in HOTSPOTS})
    engaged_labels = [HOTSPOT_LABELS[h] for h in engaged]

    lig_seq = "".join(
        protein_letters_3to1.get(r.resname, "X") for r in lig if r.id[0] == " "
    )

    return {
        "design": tag,
        "binder_len": len(lig_seq),
        "binder_seq": lig_seq,
        "n_contacts": int(n_contacts),
        "n_hbonds": int(n_hb),
        "n_saltbridges": int(n_salt),
        "dSASA_A2": round(float(dsasa), 1),
        "BSA_per_side_A2": round(float(dsasa) / 2, 1),
        "n_iface_res_receptor": len(rec_iface),
        "n_iface_res_binder": len(lig_iface),
        "hotspots_engaged": engaged_labels,
        "n_hotspots_engaged": len(engaged),
        "hotspot_coverage_pct": round(100 * len(engaged) / len(HOTSPOTS), 1),
        "iface_res_receptor": sorted(r[1] for r in rec_iface),
        "iface_res_binder": sorted(r[1] for r in lig_iface),
    }


COMPLEXES = {
    "Chai1_Best": os.path.join(DOWNLOADS, "Chai1_Best_complex.pdb"),
    "Latest_Best": os.path.join(DOWNLOADS, "Latest_Best_complex.pdb"),
    "Cyclic_Ring": os.path.join(DOWNLOADS, "Cyclic_Ring_complex.pdb"),
    "S310": os.path.join(DOWNLOADS, "S310_complex.pdb"),
    "model_4_chai": os.path.join(DOWNLOADS, "model_4_chai_complex.pdb"),
    "docked_complex": os.path.join(DOWNLOADS, "docked_complex.pdb"),
}

results = []
for tag, path in COMPLEXES.items():
    if not os.path.exists(path):
        print(f"[skip] {tag}: file not found")
        continue
    sp = split_complex(path, OUT, tag)
    if sp is None:
        print(f"[skip] {tag}: could not split")
        continue
    try:
        r = analyze(sp, tag)
    except Exception as e:
        print(f"[fail] {tag}: {e}")
        continue
    if r is None:
        print(f"[skip] {tag}: no two chains after split")
        continue
    results.append(r)
    print(f"[ok] {tag}: len={r['binder_len']} contacts={r['n_contacts']} "
          f"dSASA={r['dSASA_A2']} hotspots={r['n_hotspots_engaged']}/{len(HOTSPOTS)}")

with open(os.path.join(OUT, "interface_metrics.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"\nWrote {len(results)} interface analyses -> {OUT}")
