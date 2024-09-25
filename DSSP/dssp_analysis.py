import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.dssp import DSSP, translate
from MDAnalysis.analysis.dssp import translate, DSSP
import matplotlib.pyplot as plt
from collections import Counter

# Load your protein structure and trajectory files (replace with actual file paths)
pdb_file = 'structure.pdb'  # Replace with your PDB file
traj_file = 'eq_npt.dcd'  # Replace with your trajectory file
u = mda.Universe(pdb_file, traj_file)

# Run DSSP on the entire trajectory
dssp_analysis = DSSP(u).run()

# Get DSSP results (dssp_ndarray: 2D array, [frames x residues] of secondary structure codes)
dssp_ndarray = dssp_analysis.results.dssp_ndarray

# Translate numerical DSSP codes into readable secondary structure characters (H, E, C, etc.)
dssp_strings = [''.join(translate(dssp_frame)) for dssp_frame in dssp_ndarray]

# 1. Export Secondary Structure for Each Frame
with open("secondary_structure_per_frame.txt", "w") as f:
    for frame_idx, dssp_str in enumerate(dssp_strings):
        f.write(f"Frame {frame_idx + 1}: {dssp_str}\n")
print("Secondary structure per frame saved to 'secondary_structure_per_frame.txt'.")
