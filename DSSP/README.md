
# Protein Secondary Structure Analysis

This project analyzes the secondary structure of a protein over time using molecular dynamics simulation data. It uses DSSP (Define Secondary Structure of Proteins) algorithm to calculate the secondary structure and visualizes the results.

## Features

1. Process DSSP data from a trajectory
2. Generate a heatmap of secondary structure evolution
3. Calculate and export average secondary structure
4. Highlight specific residues in the visualization

## Dependencies

- Python 3.x
- NumPy
- Pandas
- Matplotlib
- Seaborn
- MDAnalysis

## Usage

1. Prepare your input files:
   - PDB file of your protein structure (e.g., `min.pdb`)
   - Trajectory file (e.g., `eq_npt.dcd`)

2. Run the DSSP analysis:
   ```
   python dssp_analysis.py
   ```
   This script will:
   - Perform DSSP analysis on your trajectory
   - Export secondary structure for each frame to `secondary_structure_per_frame.txt`
   - Calculate and export average secondary structure to `average_secondary_structure.txt`

3. Generate the heatmap visualization:
   ```
   python plot_heatmap.py
   ```
   This script will:
   - Read the DSSP data from `secondary_structure_per_frame.txt`
   - Create a heatmap of secondary structure evolution
   - Save the plot as `dssp_plot.png`

## File Descriptions

- `dssp_analysis.py`: Performs DSSP analysis using MDAnalysis
- `plot_heatmap.py`: Generates the heatmap visualization
- `secondary_structure_per_frame.txt`: Contains DSSP data for each frame
- `dssp_plot.png`: Heatmap visualization of secondary structure evolution

## Customization

You can customize the analysis by modifying the following:

- Input files: Update `pdb_file` and `traj_file` in `dssp_analysis.py`
- Highlighted residues: Modify `highlighted_residues` list in `plot_heatmap.py`
- Plot appearance: Adjust plot parameters in `plot_heatmap.py`

## Notes

- The heatmap uses the following color scheme:
  - 0 (blue): Coil or unstructured
  - 1 (white): β-sheet
  - 2 (red): α-helix
- Residues highlighted in red on the x-axis labels correspond to the `highlighted_residues` list

For any questions or issues, please open an issue in the repository.
