import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the fitting function
def msd_fit_function(x, D):
    return 6 * D * x

# Function to process MSD data
def process_msd(file_path, column=6):
    """
    Process MSD data to calculate diffusion coefficients for each molecule and the entire system.

    Parameters:
        file_path (str): Path to the MSD data file.
        column (int): Column to fit (2: x, 3: y, 4: z, 5: msd_axis, 6: 3D).

    Returns:
        pd.DataFrame: DataFrame with diffusion coefficients.
    """
    # Load data, skipping the first line
    data = np.loadtxt(file_path, skiprows=1)
    time = data[:, 0]  # Time is always the first column

    # Calculate number of molecules (each molecule has 5 columns, last 5 columns are system averages)
    num_molecules = (data.shape[1] - 5) // 5

    # Conversion factor from nm^2/ps to 10^-9 m^2/s
    conversion_factor = 1000

    # Initialize results
    results = []

    # Process each molecule
    for i in range(num_molecules):
        msd_data = data[:, 1 + i * 5 + (column - 2)]
        popt, _ = curve_fit(msd_fit_function, time, msd_data)
        D = popt[0] * conversion_factor  # Convert to 10^-9 m^2/s
        results.append({"Molecule": f"Molecule {i + 1}", "Diffusion Coefficient (D) [10^-9 m^2/s]": D})

        # Plot MSD and fit
        plt.figure()
        plt.plot(time, msd_data, label="MSD Data")
        plt.plot(time, msd_fit_function(time, D / conversion_factor), label=f"Fit: D = {D:.5e} 10^-9 m^2/s", linestyle="--")
        plt.xlabel("Time")
        plt.ylabel("MSD")
        plt.title(f"Molecule {i + 1} MSD")
        plt.legend()
        plt.savefig(f"msd_molecule_{i + 1}.png")
        plt.close()

    # Process system average
    system_msd_data = data[:, column]
    popt, _ = curve_fit(msd_fit_function, time, system_msd_data)
    system_D = popt[0] * conversion_factor  # Convert to 10^-9 m^2/s
    results.append({"Molecule": "System Average", "Diffusion Coefficient (D) [10^-9 m^2/s]": system_D})

    # Plot system average MSD and fit
    plt.figure()
    plt.plot(time, system_msd_data, label="System MSD Data")
    plt.plot(time, msd_fit_function(time, system_D / conversion_factor), label=f"Fit: D = {system_D:.5e} 10^-9 m^2/s", linestyle="--")
    plt.xlabel("Time")
    plt.ylabel("MSD")
    plt.title("System Average MSD")
    plt.legend()
    plt.savefig("msd_system_average.png")
    plt.close()

    # Convert results to DataFrame and save
    results_df = pd.DataFrame(results)
    results_df.to_csv("diffusion_coefficients.csv", index=False)

    return results_df

# Example usage
file_path = "msd.dat"  # Replace with the actual path to your file
column = 2  # Specify column (2: x, 3: y, 4: z, 5: msd_axis, 6: 3D)
diffusion_coefficients = process_msd(file_path, column)

# Display the results
print(diffusion_coefficients)
