# IMPORT THE PDB PARSER FROM BIOPYTHON
# This tool reads a PDB structure file and turns it into a Python object
from Bio.PDB import PDBParser

# STORE THE PATH OF OUR HIV STRUCTURE FILE IN A VARIABLE
# This makes the code easier to read and easier to change later if needed
pdb_file = 'structures/6V8X_raw.pdb'

# CREATE A PARSER OBJECT
# QUIET=True means do not print unnecessary parser warnings to the screen
parser = PDBParser(QUIET=True)

# READ THE PDB FILE AND LOAD IT INTO MEMORY AS A STRUCTURE OBJECT
# '6V8X' is just a label we give internally to this structure
structure = parser.get_structure('6V8X', pdb_file)

# PRINT A TITLE SO THE OUTPUT LOOKS CLEAN
print("Chains found in 6V8X:")

# LOOP THROUGH EACH MODEL IN THE STRUCTURE
# Most X-ray crystal structures have only one model, but this is the formal safe way
for model in structure:

    # LOOP THROUGH EACH CHAIN INSIDE THAT MODEL
    # A chain is usually labeled A, B, C, etc.
    for chain in model:

        # COUNT HOW MANY RESIDUES ARE PRESENT IN THIS CHAIN
        # Residues are amino acids like ASP, GLY, TRP etc.
        residue_count = len(list(chain))

        # PRINT THE CHAIN ID AND ITS RESIDUE COUNT
        print(f"Chain {chain.id}: {residue_count} residues")