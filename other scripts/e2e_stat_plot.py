import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def _find_e2e_column(data):
    """Helper function to find the 'End-to-End Distance' column."""
    for col in data.columns:
        if 'end-to-end distance' in col.lower():
            return col
    return None

def load_e2e_data(file_path):
    """Loads End-to-End Distance data from a single CSV file and adds a 'Molecule' column."""
    try:
        data = pd.read_csv(file_path, sep='\t')
    except pd.errors.ParserError:
        try:
            data = pd.read_csv(file_path, sep='\s+')
        except pd.errors.ParserError as e:
            raise pd.errors.ParserError(f"Could not parse CSV file: {e}")

    molecule_name = file_path.split('/')[-1].replace(".csv", "")
    data["Molecule"] = molecule_name
    return data, molecule_name

def calculate_e2e_stats(data):
    """Calculates mean and standard deviation of E2E distance."""
    e2e_column = _find_e2e_column(data)
    if e2e_column is None:
        raise KeyError("Could not find 'End-to-End Distance' column.")
    return data[e2e_column].mean(), data[e2e_column].std()

def plot_mean_e2e_with_std(stats_df, ax):
    """Plots mean E2E distance with std deviation."""
    if not all(col in stats_df.columns for col in ["Molecule", "Mean E2E", "Std E2E"]):
        raise ValueError("Stats DataFrame is missing required columns.")
    ax.bar(stats_df["Molecule"], stats_df["Mean E2E"], yerr=stats_df["Std E2E"], capsize=5, alpha=0.8)
    ax.set_title("Mean and Std Deviation of End-to-End Distance")
    ax.set_xlabel("Molecule")
    ax.set_ylabel("End-to-End Distance")
    ax.tick_params(axis='x', rotation=45)

def plot_e2e_density(combined_data, ax):
    """Plots density distribution of E2E distance."""
    e2e_column = _find_e2e_column(combined_data)
    if e2e_column is None:
        raise KeyError("Could not find 'End-to-End Distance' column in combined data.")
    sns.kdeplot(data=combined_data, x=e2e_column, hue="Molecule", ax=ax)
    ax.set_title("E2E Distance Density Plot")
    ax.set_xlabel("End-to-End Distance")
    ax.set_ylabel("Density")

def plot_e2e_time_evolution(combined_data, ax):
    """Plots time evolution of E2E distance."""
    if not all(col in combined_data.columns for col in ["Molecule", "Frame"]):
        raise ValueError("Combined DataFrame is missing required columns.")
    e2e_column = _find_e2e_column(combined_data)
    if e2e_column is None:
        raise KeyError("Could not find 'End-to-End Distance' column in combined data.")
    for molecule_name in combined_data["Molecule"].unique():
        molecule_data = combined_data[combined_data["Molecule"] == molecule_name]
        ax.plot(molecule_data["Frame"], molecule_data[e2e_column], label=molecule_name)
    ax.set_title("E2E Distance Time Evolution")
    ax.set_xlabel("Frame")
    ax.set_ylabel("End-to-End Distance")
    ax.legend()

def analyze_e2e_data(file_paths, output_dir="e2e_analysis"):
    """Analyzes End-to-End Distance data and exports plots."""
    print(f"Saving E2E analysis plots to directory: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    all_data = []
    stats = []

    for file_path in file_paths:
        data, molecule_name = load_e2e_data(file_path)
        mean_e2e, std_e2e = calculate_e2e_stats(data)
        stats.append({"Molecule": molecule_name, "Mean E2E": mean_e2e, "Std E2E": std_e2e})

    combined_data = pd.concat(all_data, ignore_index=True)
    stats_df = pd.DataFrame(stats)

    # Plot 1: Mean E2E with Std Deviation
    fig1, ax1 = plt.subplots(figsize=(8, 6))  # Adjusted size
    plot_mean_e2e_with_std(stats_df, ax1)
    fig1.tight_layout() # Add tight layout
    mean_e2e_plot_file = os.path.join(output_dir, "mean_e2e_with_std.png")
    fig1.savefig(mean_e2e_plot_file)
    plt.close(fig1)
    print(f"Mean E2E with Std Deviation plot saved to {mean_e2e_plot_file}")

    # Plot 2: E2E Density
    fig2, ax2 = plt.subplots(figsize=(8, 6))  # Adjusted size
    plot_e2e_density(combined_data, ax2)
    fig2.tight_layout() # Add tight layout
    density_plot_file = os.path.join(output_dir, "e2e_density_plot.png")
    fig2.savefig(density_plot_file)
    plt.close(fig2)
    print(f"E2E Density plot saved to {density_plot_file}")

    # Plot 3: E2E Time Evolution
    fig3, ax3 = plt.subplots(figsize=(8, 6))  # Adjusted size
    plot_e2e_time_evolution(combined_data, ax3)
    fig3.tight_layout() # Add tight layout
    time_evolution_plot_file = os.path.join(output_dir, "e2e_time_evolution.png")
    fig3.savefig(time_evolution_plot_file)
    plt.close(fig3)
    print(f"E2E Time Evolution plot saved to {time_evolution_plot_file}")

file_paths = [
    "path/to/simulation_1_e2e.csv",  # Path to the first data file
    "another/location/simulation_2_e2e.csv", # Path to the second data file
    "yet_another_e2e.csv"             # Relative path (if in the same directory)
]
analyze_e2e_data(file_paths)