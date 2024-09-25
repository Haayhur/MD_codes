import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Function to process the DSSP data file
def read_dssp_data(file_path):
    frame_data = []
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line by the colon to separate the frame number and structure
            parts = line.strip().split(': ')
            if len(parts) == 2:
                # Get the frame number (not used here, but can be useful for verification)
                frame_number = parts[0].strip().replace('Frame ', '')
                #frame_number = frame_number * 0.1
                # Get the secondary structure data (the sequence of characters)
                structure = parts[1].strip()
                # Append the structure to the list
                frame_data.append(structure)
    
    return frame_data

# Convert DSSP characters to numeric values for heatmap plotting
def dssp_to_numeric(dssp_char):
    mapping = {'E': 1, 'H': 2, '-': 0}
    return [mapping[char] for char in dssp_char]

# File path to the DSSP data
file_path = 'secondary_structure_per_frame.txt'

# Read DSSP data from the file
dssp_data = read_dssp_data(file_path)

# Convert the DSSP data to numeric form
dssp_numeric_data = [dssp_to_numeric(frame) for frame in dssp_data]

# Convert to a DataFrame for plotting
df = pd.DataFrame(dssp_numeric_data)

# Define the sequence for the x-axis
sequence = list("CPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCQYRDL")

# Assign the sequence as column labels
df.columns = sequence

# Define the residues to highlight with red color (index starts at 0)
highlighted_residues = [0, 4, 5, 8, 10, 11, 14, 18, 21, 25, 27, 34, 36, 40]

# Plot the heatmap using seaborn
plt.figure(figsize=(10, 6))
sns.heatmap(df, cmap="coolwarm", cbar_kws={'label': 'Secondary Structure'})

# Customize the x-axis labels by coloring the highlighted residues in red
ax = plt.gca()
x_labels = ax.get_xticklabels()

# Loop through the labels and change the color of the highlighted residues
for i, label in enumerate(x_labels):
    if i in highlighted_residues:
        label.set_color('red')

ax.set_xticklabels(x_labels)

plt.xlabel("Residue")
plt.ylabel("Time(ns)")
plt.title("Secondary Structure Heatmap at pH 3")
#plt.show()
plt.savefig('dssp_plot.png')
