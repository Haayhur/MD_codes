import matplotlib.pyplot as plt
import pandas as pd

def plot_rmsf_comparison(rmsf_conditions, residue_col_index=0, rmsf_col_index=1, output_image='RMSF_comparison.png'):
    """
    Plot RMSF data for multiple conditions and save the plot as an image file.

    Parameters:
    - rmsf_conditions: dict, a dictionary where keys are condition labels (e.g., 'RMSF_1')
      and values are file paths to the corresponding RMSF data CSV files.
    - residue_col_index: int, the index of the residue number column in the CSV (default is 0).
    - rmsf_col_index: int, the index of the RMSF column in the CSV (default is 1).
    - output_image: str, the filename for the saved plot image (default is 'RMSF_comparison.png').

    Output:
    - A PNG file with the RMSF comparison plot.
    """
    
    # Create a figure and axis for plotting
    fig, ax = plt.subplots(figsize=(12, 8))

    # Iterate over each condition and plot the corresponding RMSF data
    for condition, filepath in rmsf_conditions.items():
        try:
            df = pd.read_csv(filepath, sep='\t')
            ax.plot(df.iloc[:, residue_col_index], df.iloc[:, rmsf_col_index], label=condition, linewidth=2)
        except FileNotFoundError:
            print(f"File not found: {filepath}")
            continue

    # Customize the plot
    ax.set_xlabel('Residue Number', fontsize=14)
    ax.set_ylabel('RMSF (Å)', fontsize=14)
    ax.set_title('RMSF Comparison', fontsize=16, fontweight='bold')
    ax.legend(title='Condition', fontsize=12, title_fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(output_image, dpi=300, bbox_inches='tight')
    print(f"Plot saved as '{output_image}'")

def calculate_average_rmsf(rmsf_conditions, rmsf_col_index=1):
    """
    Calculate and print the average RMSF for each condition.

    Parameters:
    - rmsf_conditions: dict, a dictionary where keys are condition labels
      and values are file paths to the corresponding RMSF data CSV files.
    - rmsf_col_index: int, the index of the RMSF column in the CSV (default is 1).

    Output:
    - Prints the average RMSF for each condition.
    """
    
    print("\nAverage RMSF for each condition:")
    for condition, filepath in rmsf_conditions.items():
        try:
            df = pd.read_csv(filepath, sep='\t')
            avg_rmsf = df.iloc[:, rmsf_col_index].mean()
            print(f"{condition}: {avg_rmsf:.3f} Å")
        except FileNotFoundError:
            print(f"File not found: {filepath}")

if __name__ == "__main__":
    # Define the conditions and their corresponding file paths
    RMSF_conditions = {
        'RMSF_1': 'path/to/rmsf_1.csv',
        'RMSF_2': 'path/to/rmsf_2.csv',
        'RMSF_3': 'path/to/rmsf_3.csv'
    }

    # Plot RMSF comparison
    plot_rmsf_comparison(RMSF_conditions)

    # Calculate average RMSF for each condition
    calculate_average_rmsf(RMSF_conditions)
