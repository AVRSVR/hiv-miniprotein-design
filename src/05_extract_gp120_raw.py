# We will extract only Chain A from the 6V8X structure (raw extraction).
# "Raw" means: keep the chain A content without cleaning heteroatoms yet.

import os  # used for file/path handling
from Bio.PDB import PDBParser, PDBIO, Select  # BioPython tools to read/write PDBs

# 1) Find the downloaded 6V8X PDB/ENT file automatically in structures/
#    (Because your file name may be either 6V8X_raw.pdb or pdb6v8x.ent, etc.)
structures_dir = "structures"  # our folder that contains the PDB download(s)

# 2) Collect candidate files that look like the 6V8X input and are in PDB-like text formats
candidates = []
for fname in os.listdir(structures_dir):  # loop through all filenames inside structures/
    low = fname.lower()  # lowercase version for case-insensitive matching
    # We accept: .pdb or .ent (PDB legacy) as input
    if ("6v8x" in low) and (low.endswith(".pdb") or low.endswith(".ent")):
        candidates.append(os.path.join(structures_dir, fname))  # store full paths

# 3) Safety check: if we can't find the file, stop with a clear message
if len(candidates) == 0:
    # show what files you actually have so you can correct if needed
    existing = os.listdir(structures_dir)
    raise FileNotFoundError(
        "Could not find a 6V8X .pdb or .ent file inside structures/.\n"
        f"Files currently in structures/: {existing}\n"
    )

# 4) Choose the first candidate (usually there will be only one 6V8X raw file)
input_pdb = candidates[0]

# 5) Define the output file name for chain A extraction
output_pdb = "structures/gp120_6V8X_raw.pdb"

# 6) Define a selector that tells BioPython to KEEP ONLY chain "A"
class ChainASelect(Select):
    # accept_chain is called for every chain in the structure
    def accept_chain(self, chain):
        # keep only the gp120 reference chain (chain id == "A")
        return chain.id == "A"

# 7) Parse (read) the structure file into memory
parser = PDBParser(QUIET=True)  # suppress warnings
structure = parser.get_structure("6V8X", input_pdb)  # label = "6V8X" (internal name)

# 8) Verify chain A actually exists in the file (prevents silent mistakes)
chain_ids = set()
for model in structure:  # usually one model; loop is safe
    for chain in model:  # iterate chains
        chain_ids.add(chain.id)  # collect chain IDs

if "A" not in chain_ids:
    raise ValueError(f"Chain A not found in {input_pdb}. Found chains: {sorted(chain_ids)}")

# 9) Create an output writer and save ONLY chain A
io = PDBIO()  # BioPython object for writing structures
io.set_structure(structure)  # tell it which structure to write
io.save(output_pdb, ChainASelect())  # apply our selector and save

# 10) Print confirmation messages
print("✅ Step 5 complete: extracted Chain A (raw gp120).")
print(f"Input used: {input_pdb}")
print(f"Output saved to: {output_pdb}")