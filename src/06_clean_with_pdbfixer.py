# IMPORT PDBFIXER TO REPAIR THE STRUCTURE
from pdbfixer import PDBFixer

# IMPORT PDBFILE FROM OPENMM TO SAVE THE REPAIRED STRUCTURE
# OpenMM is the molecular physics engine that PDBFixer relies on
from openmm.app import PDBFile

# DEFINE INPUT AND OUTPUT PATHS
input_pdb = 'structures/gp120_6V8X_raw.pdb'    # The file we just extracted
output_pdb = 'structures/gp120_6V8X_clean.pdb' # The final, pristine target file

# INITIALIZE PDBFIXER WITH OUR RAW FILE
# This loads the structure into PDBFixer's memory
fixer = PDBFixer(filename=input_pdb)

# 1. FIND MISSING RESIDUES
# Looks for gaps in the amino acid sequence (e.g., jumps from residue 50 to 55)
fixer.findMissingResidues()

# 2. FIND MISSING ATOMS
# Looks for incomplete sidechains (e.g., a Lysine missing its tail)
fixer.findMissingAtoms()

# 3. ADD THE MISSING PIECES
# Uses physics-based modeling to rebuild the missing residues and atoms
fixer.addMissingAtoms()

# 4. REMOVE HETEROATOMS
# keepWater=False ensures we strip out any lingering water molecules or ions
# We only want pure protein for the design phase
fixer.removeHeterogens(keepWater=False)

# 5. ADD MISSING HYDROGENS AT pH 7.4
# CRITICAL STEP: Calculates the correct protonation states for all amino acids
# at human physiological pH. This allows docking software to calculate charges.
fixer.addMissingHydrogens(7.4)

# 6. SAVE THE CLEANED STRUCTURE
# Opens a new file in write mode ('w')
with open(output_pdb, 'w') as f:
    # Writes the repaired topology (atoms/bonds) and positions (3D coordinates)
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

# PRINT CONFIRMATION
print("✅ Step 6 complete: Structure repaired and cleaned.")
print(f"💾 Saved pristine target to: {output_pdb}")
print("🧪 Hydrogens added at pH 7.4. Ready for hotspot mapping.")