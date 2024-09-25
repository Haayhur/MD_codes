import pandas as pd
import matplotlib.pyplot as plt

def plot_data(file_path, label, frame_col="Frame", rg_col="Radius of Gyration (A)", conversion_factor=10e-6):
    """
    Plot radius of gyration data over time from a CSV file.

    Parameters:
    - file_path: str, path to the CSV file containing the data.
    - label: str, label for the plot line.
    - frame_col: str, column name for the frame numbers in the CSV (default is 'Frame').
    - rg_col: str, column name for the radius of gyration values in the CSV (default is 'Radius of Gyration (A)').
    - conversion_factor: float, conversion factor from frames to time units (default is 10e-6 ms).
    
    Output:
    - Plots the radius of gyration data against time.
    """
    
    # Read the data into a pandas DataFrame
    data = pd.read_csv(file_path, header=None, names=[frame_col, rg_col])
    
    # Convert columns to numeric, coercing errors into NaN (e.g., for non-numeric values)
    data[frame_col] = pd.to_numeric(data[frame_col], errors='coerce')
    data[rg_col] = pd.to_numeric(data[rg_col], errors='coerce')
    
    # Convert frames to time (e.g., microseconds)
    data["Time (ms)"] = data[frame_col] * conversion_factor
    
    # Plot the radius of gyration over time
    plt.plot(data["Time (ms)"], data[rg_col], label=label)

if __name__ == "__main__":
    # Define file paths and labels for the different Rg data
    file_paths = {
        'Rg_1': 'path/to/ligand_rg_values_1.csv',
        'Rg_2': 'path/to/ligand_rg_values_2.csv',
        'Rg_3': 'path/to/ligand_rg_values_3.csv'
    }

    # Iterate through each file and plot its data
    for label, file_path in file_paths.items():
        try:
            plot_data(file_path, label)
        except FileNotFoundError:
            print(f"File not found: {file_path}")

    # Customize the plot
    plt.xlabel('Time (ms)', fontsize=12)
    plt.ylabel('Radius of Gyration (Å)', fontsize=12)
    plt.title('Radius of Gyration Over Time', fontsize=14, fontweight='bold')
    plt.legend(title='Condition', fontsize=10, title_fontsize=12)
    plt.grid(True)

    # Save the plot to a file
    plt.savefig('Rg_comparison.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'Rg_comparison.png'")
