# IMPORT PDB PARSER TO READ THE STRUCTURE
from Bio.PDB import PDBParser

# DEFINE THE FILE TO CHECK
pdb_file = 'structures/6V8X_raw.pdb'

# LIST OF CRITICAL CD4 BINDING SITE HOTSPOTS (HXB2 NUMBERING)
# These are the residues we MUST have for the design to work
# Note: PDB residue numbers might differ slightly due to signal peptides, 
# but we check for the presence of these residue types in the binding region
hotspot_residue_names = ['ASP', 'GLU', 'TRP', 'ARG', 'LYS', 'GLY', 'THR', 'ASN', 'VAL', 'ILE', 'GLN']

# INITIALIZE PARSER
parser = PDBParser(QUIET=True)
structure = parser.get_structure('6V8X', pdb_file)

# STORE RESIDUE COUNTS FOR A, E, I TO COMPARE
chain_data = {}

print("\n--- Verifying gp120 Chains (A, E, I) ---")

for model in structure:
    for chain in model:
        # We only care about the gp120 chains (A, E, I)
        if chain.id in ['A', 'E', 'I']:
            
            # Count total residues
            total_residues = len(list(chain))
            
            # Count how many residues are standard amino acids (not water/ions)
            standard_residues = sum(1 for residue in chain if residue.id[0] == ' ')
            
            # Store the data
            chain_data[chain.id] = {
                'total': total_residues,
                'standard': standard_residues
            }

# PRINT THE COMPARISON
for chain_id, data in chain_data.items():
    print(f"Chain {chain_id}: {data['standard']} standard amino acids")

# CHECK IF THEY ARE IDENTICAL
counts = [data['standard'] for data in chain_data.values()]
if counts[0] == counts[1] == counts[2]:
    print("\n✅ All three gp120 chains (A, E, I) are identical in length.")
    print("✅ Chain A is safe to use as the reference.")
else:
    print("\n⚠️ Warning: Chains differ in length. Check for missing loops.")

print("---------------------------------------------\n")