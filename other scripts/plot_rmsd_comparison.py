import matplotlib.pyplot as plt
import pandas as pd

def plot_rmsd_comparison(rmsd_conditions, time_col_index=0, rmsd_col_index=1, output_image='RMSD_comparison.png'):
    """
    Plot RMSD data for multiple conditions and save the plot as an image file.

    Parameters:
    - rmsd_conditions: dict, a dictionary where keys are condition labels (e.g., 'RMSD_1')
      and values are file paths to the corresponding RMSD data CSV files.
    - time_col_index: int, the index of the time column in the CSV (default is 0).
    - rmsd_col_index: int, the index of the RMSD column in the CSV (default is 1).
    - output_image: str, the filename for the saved plot image (default is 'RMSD_comparison.png').

    Output:
    - A PNG file with the RMSD comparison plot.
    """
    
    # Create a figure and axis for plotting
    fig, ax = plt.subplots(figsize=(12, 8))

    # Iterate over each condition and plot the corresponding RMSD data
    for condition, filepath in rmsd_conditions.items():
        try:
            df = pd.read_csv(filepath)
            ax.plot(df.iloc[:, time_col_index], df.iloc[:, rmsd_col_index], label=condition, linewidth=2)
        except FileNotFoundError:
            print(f"File not found: {filepath}")
            continue

    # Customize the plot
    ax.set_xlabel('Time (ps)', fontsize=14)
    ax.set_ylabel('RMSD (Å)', fontsize=14)
    ax.set_title('RMSD Comparison', fontsize=16, fontweight='bold')
    ax.legend(title='Condition', fontsize=12, title_fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"Plot saved as '{output_image}'")

def calculate_average_rmsd(rmsd_conditions, rmsd_col_index=1):
    """
    Calculate and print the average RMSD for each condition.

    Parameters:
    - rmsd_conditions: dict, a dictionary where keys are condition labels
      and values are file paths to the corresponding RMSD data CSV files.
    - rmsd_col_index: int, the index of the RMSD column in the CSV (default is 1).

    Output:
    - Prints the average RMSD for each condition.
    """
    
    print("\nAverage RMSD for each condition:")
    for condition, filepath in rmsd_conditions.items():
        try:
            df = pd.read_csv(filepath)
            avg_rmsd = df.iloc[:, rmsd_col_index].mean()
            print(f"{condition}: {avg_rmsd:.3f} Å")
        except FileNotFoundError:
            print(f"File not found: {filepath}")

if __name__ == "__main__":
    # Define the conditions and their corresponding file paths
    RMSD_conditions = {
        'RMSD_1': 'path/to/rmsd_1.csv',
        'RMSD_2': 'path/to/rmsd_2.csv',
        'RMSD_3': 'path/to/rmsd_3.csv'
    }

    # Plot RMSD comparison
    plot_rmsd_comparison(RMSD_conditions)

    # Calculate average RMSD for each condition
    calculate_average_rmsd(RMSD_conditions)
