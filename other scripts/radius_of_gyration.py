import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def analyze_radius_of_gyration(pdb_file, dcd_file, selection="resnum 1:225"):
    """
    Analyze the radius of gyration (Rg) for a selected group of atoms
    from a molecular dynamics trajectory and save the results as a plot and CSV file.

    Parameters:
    - pdb_file: str, the path to the PDB file (e.g., 'min.pdb')
    - dcd_file: str, the path to the DCD trajectory file (e.g., 'eq_npt.dcd')
    - selection: str, the atom selection string (default is 'resnum 1:225')

    Output:
    - A CSV file containing frame numbers and Rg values.
    - A PNG file with the Rg plot over time.
    """

    # Initialize the MDAnalysis universe
    u = mda.Universe(pdb_file, dcd_file)

    # Select the atoms based on the provided selection (default: resnum 1:225)
    ligand = u.select_atoms(selection)

    # Lists to store frame numbers and Rg values
    frames = []
    ligand_rg_values = []

    # Iterate through all frames and calculate Rg
    for ts in u.trajectory:
        ligand_rgyr = ligand.radius_of_gyration()
        frames.append(ts.frame)
        ligand_rg_values.append(ligand_rgyr)

    # Convert to numpy arrays for analysis
    frames = np.array(frames)
    ligand_rg_values = np.array(ligand_rg_values)

    # Calculate statistics
    avg_rg = np.mean(ligand_rg_values)
    std_rg = np.std(ligand_rg_values)

    # Print statistics
    print(f"\nLigand Statistics:\nAverage Rg: {avg_rg:.2f} Å\nStandard Deviation: {std_rg:.2f} Å")

    # Create a DataFrame and export to CSV
    df = pd.DataFrame({"Frame": frames, "Radius of Gyration (Å)": ligand_rg_values})
    df.to_csv("ligand_rg_values.csv", index=False)

    print("Data has been exported to 'ligand_rg_values.csv'.")

    return frames, ligand_rg_values, avg_rg, std_rg

def plot_radius_of_gyration(frames, rg_values):
    """
    Plot the radius of gyration (Rg) over time.

    Parameters:
    - frames: array-like, the frame numbers
    - rg_values: array-like, the corresponding Rg values
    """
    plt.figure(figsize=(10, 6))
    plt.plot(frames, rg_values, label='Ligand', color='b')
    plt.xlabel('Frame')
    plt.ylabel('Radius of Gyration (Å)')
    plt.title('Radius of Gyration over Time')
    plt.legend()
    plt.tight_layout()
    plt.savefig('rg_over_time.png')
    print("Plot saved as 'rg_over_time.png'.")

if __name__ == "__main__":
    pdb_file = "min.pdb"
    dcd_file = "eq_npt.dcd"
    
    # Run the analysis
    frames, ligand_rg_values, avg_rg, std_rg = analyze_radius_of_gyration(pdb_file, dcd_file)
    
    # Plot the radius of gyration
    plot_radius_of_gyration(frames, ligand_rg_values)
