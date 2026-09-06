# IMPORTING THE BIOLOGICAL DATABASE TOOL FROM BIOPYTHON LIBRARY
# BioPython is a collection of tools for analyzing protein sequences and structures
from Bio.PDB import PDBList
import os

# INITIALIZES THE PDB LIST CONNECTION OBJECT
# This object knows how to download files from the RCSB Protein Data Bank
pdbl = PDBList()

# DOWNLOADS THE ACTUAL PDB FILE FROM THE INTERNET
# '6V8X' is the specific structure ID from the Protein Data Bank
# pdir='./structures' tells it to put the file inside our structures folder
# overwrite=True ensures we get the latest version if we run this twice
# file_format='pdb' forces download as .pdb (not .cif) - most of our tools expect .pdb
pdbl.retrieve_pdb_file('6V8X', pdir='./structures', overwrite=True, file_format='pdb')

# PRINTS A CONFIRMATION MESSAGE TO THE TERMINAL SO WE KNOW IT WORKED
print("✅ Downloaded 6V8X successfully")

# PRINTS THE PATH WHERE THE FILE WAS SAVED ON YOUR COMPUTER
print("📁 Saved to: ./structures/")