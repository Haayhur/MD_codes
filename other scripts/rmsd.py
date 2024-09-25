import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt
import numpy as np

def calculate_rmsd(trajectory_file, topology_file, selection="resnum 1:225"):
    """
    Calculate the RMSD of selected atoms from a molecular dynamics trajectory
    and save the results as a plot and CSV file.

    Parameters:
    - trajectory_file: str, the path to the trajectory file (e.g., 'eq_npt.dcd')
    - topology_file: str, the path to the topology file (e.g., 'min.pdb')
    - selection: str, the atom selection string for the ligand or region of interest (default is 'resnum 1:225')

    Output:
    - A PNG file showing the RMSD plot over time.
    - A CSV file containing the time and RMSD data.
    """

    # Load the trajectory and topology
    universe = mda.Universe(topology_file, trajectory_file)

    # Select the atoms based on the provided selection
    ligand = universe.select_atoms(selection)

    # Create a reference from the first frame and align it to the center of mass
    reference = ligand.positions.copy()
    ligand.translate(-ligand.center_of_mass())

    # Create the RMSD analyzer
    rmsd_analyzer = rms.RMSD(ligand, reference=reference)
    
    # Run the RMSD calculation
    rmsd_analyzer.run()

    # Get the results
    time = rmsd_analyzer.times
    rmsd = rmsd_analyzer.rmsd[:, 2]  # The third column contains the RMSD values

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(time, rmsd, label="RMSD")
    plt.xlabel("Time (ps)")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD over Time")
    plt.legend()
    plt.tight_layout()
    plt.savefig("rmsd.png")
    
    # Save RMSD data to CSV
    np.savetxt("rmsd_data.csv", np.column_stack((time, rmsd)), delimiter=",", 
               header="Time (ps),RMSD (A)", comments="")
    
    # Print average RMSD
    avg_rmsd = np.mean(rmsd)
    print(f"Average RMSD: {avg_rmsd:.2f} Å")

    return avg_rmsd

# Run the function
average_rmsd = calculate_rmsd("eq_npt.dcd", "min.pdb")
