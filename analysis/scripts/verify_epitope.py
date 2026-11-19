"""
Independent check of WHERE each design binds on gp120.

For each complex, compute the CD4-binding-site centroid from that complex's own
hotspot residues (HXB2 numbering), then measure the distance from the binder's
center of mass to that centroid, and the minimum binder-to-hotspot distance.
This is independent of the contact-based hotspot call in interface_analysis.py.
"""
import os, glob, json, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SPLIT_DIR = os.path.join(PROJECT, "paper_analysis", "01_structures")

# core CD4-binding-site residues (HXB2): CD4 binding loop, bridging sheet, V5 contact
CD4BS_CORE = [368, 370, 371, 425, 427, 429, 457, 473, 280]

parser = PDBParser(QUIET=True)
rows = []

for f in sorted(glob.glob(os.path.join(SPLIT_DIR, "*_split.pdb"))):
    tag = os.path.basename(f).replace("_split.pdb", "")
    st = parser.get_structure(tag, f)
    m = next(iter(st))
    if "A" not in m or "B" not in m:
        continue
    rec, lig = m["A"], m["B"]

    hs_coords = []
    found = []
    for res in rec:
        if res.id[0] != " ":
            continue
        if res.id[1] in CD4BS_CORE and "CA" in res:
            hs_coords.append(res["CA"].coord)
            found.append(res.id[1])
    if len(hs_coords) < 3:
        print(f"[skip] {tag}: only {len(hs_coords)} CD4bs residues found")
        continue
    hs_coords = np.array(hs_coords)
    cd4bs_centroid = hs_coords.mean(axis=0)

    lig_atoms = np.array([a.coord for a in lig.get_atoms() if a.element != "H"])
    lig_com = lig_atoms.mean(axis=0)

    dist_com = float(np.linalg.norm(lig_com - cd4bs_centroid))
    # closest approach of any binder atom to any CD4bs hotspot CA
    dmat = np.linalg.norm(lig_atoms[:, None, :] - hs_coords[None, :, :], axis=-1)
    min_dist = float(dmat.min())
    # how many hotspots have a binder atom within 8 A of their CA
    near = int((dmat.min(axis=0) < 8.0).sum())

    rows.append({
        "design": tag,
        "cd4bs_residues_found": len(found),
        "binder_COM_to_CD4bs_centroid_A": round(dist_com, 1),
        "min_binder_to_hotspot_A": round(min_dist, 1),
        "hotspots_within_8A": near,
        "n_hotspots_checked": len(found),
        "verdict": "CD4bs-engaged" if near >= 4 else ("peripheral" if near >= 1 else "off-site"),
    })

rows.sort(key=lambda r: r["binder_COM_to_CD4bs_centroid_A"])
print(f"{'DESIGN':<18}{'COM->CD4bs(A)':>15}{'min dist(A)':>13}{'HS<8A':>8}  VERDICT")
print("-" * 70)
for r in rows:
    print(f"{r['design']:<18}{r['binder_COM_to_CD4bs_centroid_A']:>15}"
          f"{r['min_binder_to_hotspot_A']:>13}{r['hotspots_within_8A']:>4}/{r['n_hotspots_checked']:<3}  {r['verdict']}")

with open(os.path.join(SPLIT_DIR, "epitope_verification.json"), "w") as f:
    json.dump(rows, f, indent=2)
