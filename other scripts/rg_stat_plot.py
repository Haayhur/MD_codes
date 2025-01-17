import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_rg_data(file_path):
    """Loads Rg data from a single CSV file and adds a 'Molecule' column."""
    try:
        data = pd.read_csv(file_path, sep='\t')  # Assuming tab-separated, adjust if needed
    except pd.errors.ParserError:
        try:
            data = pd.read_csv(file_path, sep='\s+')  # Trying whitespace as separator
        except pd.errors.ParserError as e:
            raise pd.errors.ParserError(f"Could not parse CSV file with tab or whitespace delimiters: {e}")

    molecule_name = file_path.split('/')[-1].replace(".csv", "")
    data["Molecule"] = molecule_name
    return data, molecule_name

def calculate_rg_stats(data):
    """Calculates mean and standard deviation of Rg for a given DataFrame."""
    # Determine the correct column name for Radius of Gyration
    rg_column = None
    for col in data.columns:
        if 'Radius of Gyration' in col:
            rg_column = col
            break
    if rg_column is None:
        raise KeyError("Could not find 'Radius of Gyration' column in the CSV file.")

    return data[rg_column].mean(), data[rg_column].std()

def plot_mean_rg_with_std(stats_df, ax):
    """Plots the mean Rg with standard deviation for each molecule."""
    ax.bar(stats_df["Molecule"], stats_df["Mean Rg"], yerr=stats_df["Std Rg"], capsize=5, alpha=0.8)
    ax.set_title("Mean and Std Deviation of Rg")
    ax.set_xlabel("Molecule")
    ax.set_ylabel("Radius of Gyration")  # Removed (A) for generality
    ax.tick_params(axis='x', rotation=45)

def plot_rg_density(combined_data, ax):
    """Plots the density distribution of Rg for each molecule."""
    # Determine the correct column name for Radius of Gyration
    rg_column = None
    for col in combined_data.columns:
        if 'Radius of Gyration' in col:
            rg_column = col
            break
    if rg_column is None:
        raise KeyError("Could not find 'Radius of Gyration' column in the combined data.")

    sns.kdeplot(data=combined_data, x=rg_column, hue="Molecule", ax=ax)
    ax.set_title("Rg Density Plot")
    ax.set_xlabel("Radius of Gyration") # Removed (A) for generality
    ax.set_ylabel("Density")

def plot_rg_time_evolution(combined_data, ax):
    """Plots the time evolution of Rg for each molecule."""
    # Determine the correct column name for Radius of Gyration
    rg_column = None
    for col in combined_data.columns:
        if 'Radius of Gyration' in col:
            rg_column = col
            break
    if rg_column is None:
        raise KeyError("Could not find 'Radius of Gyration' column in the combined data.")

    for molecule_name in combined_data["Molecule"].unique():
        molecule_data = combined_data[combined_data["Molecule"] == molecule_name]
        ax.plot(molecule_data["Frame"], molecule_data[rg_column], label=molecule_name)
    ax.set_title("Rg Time Evolution")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Radius of Gyration") # Removed (A) for generality
    ax.legend()


def analyze_rg_data(file_paths, output_dir="rg_analysis"):
    """
    Analyzes Radius of Gyration (Rg) data from multiple simulations and exports each plot separately.

    Parameters:
        file_paths (list of str): List of CSV file paths containing Rg data.
        output_dir (str): Directory to save the individual plots.

    Each CSV file should have the following columns:
        - Frame: Frame number.
        - Time (ps): Time in picoseconds.
        - Radius of Gyration: Rg values for the corresponding frame.
    """
    import os

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    all_data = []
    stats = []

    for file_path in file_paths:
        data, molecule_name = load_rg_data(file_path)
        mean_rg, std_rg = calculate_rg_stats(data)
        stats.append({"Molecule": molecule_name, "Mean Rg": mean_rg, "Std Rg": std_rg})
        all_data.append(data)

    combined_data = pd.concat(all_data, ignore_index=True)
    stats_df = pd.DataFrame(stats)

    # Plot 1: Mean Rg with Std Deviation
    fig1, ax1 = plt.subplots(figsize=(7, 5))
    plot_mean_rg_with_std(stats_df, ax1)
    mean_rg_plot_file = os.path.join(output_dir, "mean_rg_with_std.png")
    fig1.savefig(mean_rg_plot_file)
    plt.close(fig1)
    print(f"Mean Rg with Std Deviation plot saved to {mean_rg_plot_file}")

    # Plot 2: Rg Density
    fig2, ax2 = plt.subplots(figsize=(7, 5))
    plot_rg_density(combined_data, ax2)
    density_plot_file = os.path.join(output_dir, "rg_density_plot.png")
    fig2.savefig(density_plot_file)
    plt.close(fig2)
    print(f"Rg Density plot saved to {density_plot_file}")

    # Plot 3: Rg Time Evolution
    fig3, ax3 = plt.subplots(figsize=(7, 5))
    plot_rg_time_evolution(combined_data, ax3)
    time_evolution_plot_file = os.path.join(output_dir, "rg_time_evolution.png")
    fig3.savefig(time_evolution_plot_file)
    plt.close(fig3)
    print(f"Rg Time Evolution plot saved to {time_evolution_plot_file}")

# Example usage
file_paths = [
    "fda1_rog.csv",
    "fda_rog.csv",
]
analyze_rg_data(file_paths)