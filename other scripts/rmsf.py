import MDAnalysis as mda
from MDAnalysis.analysis import rms
import pandas as pd
import matplotlib.pyplot as plt

def calculate_rmsf(trajectory_file, topology_file, selection="protein"):
    """
    Calculate the RMSF (Root Mean Square Fluctuation) for selected atoms
    in a molecular dynamics trajectory and save the results as a plot and CSV file.

    Parameters:
    - trajectory_file: str, the path to the trajectory file (e.g., 'eq_npt.dcd')
    - topology_file: str, the path to the topology file (e.g., 'min.pdb')
    - selection: str, the atom selection string (default is 'protein')

    Output:
    - A PNG file with the RMSF plot.
    - A CSV file containing residue numbers and their corresponding RMSF values.
    """

    # Load the trajectory and topology
    universe = mda.Universe(topology_file, trajectory_file)

    # Select atoms for RMSF calculation (default: alpha carbons of the protein)
    atoms = universe.select_atoms(selection)

    # Run RMSF calculation
    rmsf_analysis = rms.RMSF(atoms, verbose=True).run()

    # Create a DataFrame with residue numbers and RMSF values
    rmsf_data = pd.DataFrame({'Residue': atoms.resnums, 'RMSF': rmsf_analysis.results.rmsf})

    # Save RMSF data to CSV
    rmsf_data.to_csv("rmsf.csv", sep='\t', index=False, encoding='utf8')

    # Plot RMSF
    plt.figure(figsize=(12, 6))
    plt.plot(rmsf_data['Residue'], rmsf_data['RMSF'], label='RMSF')
    plt.xlabel('Residue Number')
    plt.ylabel('RMSF (Å)')
    plt.title('RMSF of Protein')
    plt.legend()
    plt.tight_layout()
    plt.savefig('RMSF_plot.png')

    # Print average RMSF
    avg_rmsf = rmsf_data['RMSF'].mean()
    print(f"Average RMSF: {avg_rmsf:.2f} Å")

    return avg_rmsf

# Run the function
average_rmsf = calculate_rmsf("eq_npt.dcd", "min.pdb")
