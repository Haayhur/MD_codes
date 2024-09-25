import numpy as np
import matplotlib.pyplot as plt

# Function to calculate the degree of deprotonation based on pH, pKa, and charge q
def calculate_degree_of_deprotonation(ph_values, pka, q=1):
    """
    Calculate the degree of deprotonation of a compound at different pH values.
    
    Parameters:
    ph_values (numpy array): Array of pH values.
    pka (float): The pKa value of the compound.
    q (int, optional): Charge value. Default is 1.
    
    Returns:
    numpy array: Degree of deprotonation for each pH value.
    """
    return 1 / (10**(q * (pka - ph_values)) + 1)

# Example pKa of the compound (can be adjusted)
pka = 6.5
# Charge value (default is 1 for simplicity)
q = 1

# Generate pH values from 3 to 12 (can be adjusted based on the compound's range)
ph_values = np.arange(3, 12, 0.1)

# Calculate degree of deprotonation for each pH value
degree_of_deprotonation = calculate_degree_of_deprotonation(ph_values, pka, q)

# Calculate degree of protonation (1 - degree of deprotonation)
degree_of_protonation = 1 - degree_of_deprotonation

# Specific pH values to highlight
highlight_pH = [3, 5.5, 7.4]
highlight_degree = calculate_degree_of_deprotonation(np.array(highlight_pH), pka, q)

# Plot the results
plt.plot(ph_values, degree_of_protonation, label='Degree of Protonation')
plt.scatter(highlight_pH, 1 - highlight_degree, color='red', zorder=5, label='pH of Interest')

# Annotate the highlighted points
for i, ph in enumerate(highlight_pH):
    plt.text(ph, 1 - highlight_degree[i], f'pH {ph}', color='black', ha='right')

plt.title('Degree of Protonation vs pH for a Compound')
plt.xlabel('pH')
plt.ylabel('Degree of Protonation')
plt.grid(True)
plt.legend()

# Save the plot as a PNG file
plt.savefig("pH_titration.png")
