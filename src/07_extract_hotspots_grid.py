import numpy as np
import json
from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner

# ============================================================
# INPUT / OUTPUT
# ============================================================
input_pdb = 'structures/gp120_6V8X_clean.pdb'
output_json = 'structures/cd4bs_hotspots_6V8X.json'

# ============================================================
# HXB2 GP120 SEQUENCE ONLY (residues 31-511, 481 residues)
# This is the EXACT gp120 region from HXB2 (UniProt P04578)
# ============================================================
hxb2_gp120 = (
    "TEKLWVTVYYGVPVWKDAETTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLENVTENFNMWKNNMVEQMHEDIISLWDQ"
    "SLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCN"
    "TSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFT"
    "DNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKIREQFGNNTII"
    "IFNQSSGGDPEIVMHSFNCGGEFFYCNTSQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQEVGKAMYAPP"
    "IRGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"
)

print(f"HXB2 gp120 length: {len(hxb2_gp120)}")

# ============================================================
# HOTSPOTS (relative to HXB2 gp120 start at residue 31)
# So HXB2 368 = position (368-31) = 337 in this string
# ============================================================
hotspots = [
    {"label": "D368", "pos": 368, "expected": "D"},  # 368-31=337
    {"label": "E370", "pos": 370, "expected": "E"},  # 370-31=339
    {"label": "W427", "pos": 427, "expected": "W"},  # 427-31=396
    {"label": "D457", "pos": 457, "expected": "D"},  # 457-31=426
    {"label": "R429", "pos": 429, "expected": "R"},  # 429-31=398
    {"label": "K421", "pos": 421, "expected": "K"},  # 421-31=390
    {"label": "G473", "pos": 473, "expected": "G"},  # 473-31=442
    {"label": "E461", "pos": 461, "expected": "E"},  # 461-31=430
    {"label": "T257", "pos": 257, "expected": "T"},  # 257-31=226
    {"label": "N280", "pos": 280, "expected": "N"},  # 280-31=249
    {"label": "V255", "pos": 255, "expected": "V"},  # 255-31=224
    {"label": "I371", "pos": 371, "expected": "I"},  # 371-31=340
    {"label": "Q105", "pos": 105, "expected": "Q"},  # 105-31=74
    {"label": "N425", "pos": 425, "expected": "N"},  # 425-31=394
]

# ============================================================
# HELPERS
# ============================================================
three_to_one = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}

# ============================================================
# 1. EXTRACT PDB SEQUENCE
# ============================================================
parser = PDBParser(QUIET=True)
structure = parser.get_structure('gp120', input_pdb)
model = next(structure.get_models())
chain = model['A']

pdb_residues = []
pdb_seq = []
for r in chain:
    if r.id[0] == ' ':
        name3 = r.get_resname()
        name1 = three_to_one.get(name3, 'X')
        pdb_residues.append({'res_obj': r, 'name3': name3, 'name1': name1, 'pdb_num': r.id[1]})
        pdb_seq.append(name1)

pdb_seq_str = ''.join(pdb_seq)
print(f"PDB sequence length: {len(pdb_seq_str)}")

# ============================================================
# 2. GLOBAL ALIGNMENT (HXB2 gp120 vs PDB sequence)
# ============================================================
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -1

best = aligner.align(hxb2_gp120, pdb_seq_str)[0]
print(f"Alignment score: {best.score}")

# ============================================================
# 3. BUILD MAPPING: HXB2 position -> PDB residue index
# ============================================================
# Convert HXB2 hotspot positions to indices in the HXB2 gp120 string
# HXB2 gp120 starts at residue 31, so subtract 31 to get string index
hxb2_to_pdb_idx = {}  # maps HXB2 absolute position -> index in pdb_residues list

# Walk through the alignment block by block
target_blocks = best.aligned[0]  # HXB2 gp120 blocks
query_blocks = best.aligned[1]   # PDB sequence blocks

hxb2_pos = 31  # Start of HXB2 gp120
pdb_idx = 0

for block_idx in range(len(target_blocks)):
    t_start, t_end = target_blocks[block_idx]
    q_start, q_end = query_blocks[block_idx]
    
    # Within this block, each aligned position maps 1-to-1
    block_len = t_end - t_start
    for i in range(block_len):
        abs_hxb2_pos = hxb2_pos + t_start + i - 31  # Wait, no
        # Actually, t_start is the index in hxb2_gp120 string
        # The absolute HXB2 position is 31 + t_start + i
        # The PDB index is q_start + i
        hxb2_abs_pos = 31 + t_start + i
        hxb2_to_pdb_idx[hxb2_abs_pos] = q_start + i
    
    hxb2_pos += block_len

print(f"Mapped {len(hxb2_to_pdb_idx)} HXB2 positions to PDB indices")

# ============================================================
# 4. EXTRACT COORDINATES FOR ALL 14 HOTSPOTS
# ============================================================
hotspot_coords = []
hotspot_records = []

print("\n--- Hotspot mapping results ---")
for h in hotspots:
    hxpos = h['pos']
    label = h['label']
    expected = h['expected']
    
    if hxpos not in hxb2_to_pdb_idx:
        print(f"❌ {label}: HXB2 pos {hxpos} NOT MAPPED")
        continue
    
    pdb_idx = hxb2_to_pdb_idx[hxpos]
    
    if pdb_idx >= len(pdb_residues):
        print(f"❌ {label}: PDB index {pdb_idx} out of range")
        continue
    
    entry = pdb_residues[pdb_idx]
    res_obj = entry['res_obj']
    name3 = entry['name3']
    name1 = entry['name1']
    pdb_num = entry['pdb_num']
    
    type_match = (name1 == expected)
    
    try:
        ca = res_obj['CA']
        coord = ca.get_coord()
        hotspot_coords.append(coord)
        
        hotspot_records.append({
            "hotspot": label,
            "hxb2_pos": hxpos,
            "expected": expected,
            "found": name1,
            "found_3letter": name3,
            "match": type_match,
            "pdb_resseq": pdb_num,
            "ca_xyz": [float(coord[0]), float(coord[1]), float(coord[2])]
        })
        
        status = "✓" if type_match else "⚠ TYPE DIFFERS"
        print(f"{status} {label}: HXB2 {hxpos} -> PDB idx {pdb_idx}, resseq {pdb_num}, found {name3} ({name1})")
    except KeyError:
        print(f"❌ {label}: Missing CA atom at PDB idx {pdb_idx}")

# ============================================================
# 5. CENTROID AND POCKET DIAMETER
# ============================================================
if len(hotspot_coords) < 14:
    print(f"\n❌ Only found {len(hotspot_coords)}/14 hotspots. Cannot compute grid.")
    print("This likely means the alignment didn't cover all hotspot positions.")
    print("Check if your PDB sequence covers the full CD4bs region.")
    exit(1)

coords_arr = np.array(hotspot_coords)
centroid = coords_arr.mean(axis=0)
dist_mat = np.linalg.norm(coords_arr[:, None, :] - coords_arr[None, :, :], axis=-1)
pocket_diam = float(dist_mat.max())
box_size = pocket_diam + 20.0

# ============================================================
# 6. SAVE JSON
# ============================================================
grid_data = {
    "structure": "6V8X",
    "chain": "A",
    "method": "global alignment of HXB2 gp120 (residues 31-511) to PDB sequence",
    "centroid": {"x": float(centroid[0]), "y": float(centroid[1]), "z": float(centroid[2])},
    "pocket_diameter_A": pocket_diam,
    "box_size_A": box_size,
    "hotspots": hotspot_records
}

with open(output_json, 'w') as f:
    json.dump(grid_data, f, indent=2)

print(f"\n{'='*60}")
print("✅ STEP 7 COMPLETE")
print(f"   {len(hotspot_coords)}/14 hotspots mapped successfully")
print(f"   Centroid: {centroid[0]:.3f}  {centroid[1]:.3f}  {centroid[2]:.3f}")
print(f"   Pocket diameter: {pocket_diam:.2f} Å")
print(f"   Box size: {box_size:.2f} Å")
print(f"   Saved to: {output_json}")
print(f"{'='*60}")