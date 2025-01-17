import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_diffusion_coefficients(input_files, highlights, output_file="diffusion_coefficients_plot.png"):
    """
    Plot the diffusion coefficients for each molecule from multiple input MSD data files using a line plot.

    Parameters:
        input_files (list of str): List of paths to the CSV files containing diffusion coefficients.
        highlights (dict): Dictionary where keys are filenames and values are lists of molecule indices to highlight.
        output_file (str): Path to save the output plot.
    """
    # Initialize lists for data
    all_molecules = []
    all_diffusion_coefficients = []

    # Process each file
    for file in input_files:
        data = pd.read_csv(file)
        molecules = data["Molecule"]
        diffusion_coefficients = data["Diffusion Coefficient (D) [10^-9 m^2/s]"]

        all_molecules.extend([f"{mol} ({file})" for mol in molecules])
        all_diffusion_coefficients.append(diffusion_coefficients)

    # Generate the plot
    plt.figure(figsize=(12, 8))

    # Plot each file's data
    for i, file in enumerate(input_files):
        coefficients = all_diffusion_coefficients[i]
        indices = range(len(coefficients))

        # Plot all points in blue by default
        plt.plot(indices, coefficients, marker='o', label=f"File: {file}", linestyle='-', color='blue')

        # Highlight specified points in red
        if file in highlights:
            highlight_indices = highlights[file]
            plt.scatter(highlight_indices, [coefficients[j] for j in highlight_indices], color='red', label=f"Highlighted in {file}", zorder=5)

    # Add labels and title
    plt.xlabel("Molecule Index")
    plt.ylabel("Diffusion Coefficient (D) [10^-9 m²/s]")
    plt.title("Diffusion Coefficients for Each Molecule")

    # Add legend
    plt.legend()

    # Save the plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    print(f"Plot saved to {output_file}")

# Example usage
# Replace ['file1.csv', 'file2.csv'] with the paths to your CSV files
# Specify highlights as {filename: [indices]}
plot_diffusion_coefficients(
    ["diffusion_coefficients.csv"],
    highlights={"diffusion_coefficients.csv": [1, 7]}
)
