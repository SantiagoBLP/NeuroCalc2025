#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def plot_bands():
    """Plot band structure from Quantum Espresso output and save to a file."""
    try:
        # Determine the directory where the script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Load data from file 'bands.dat.gnu' located in the same directory
        data_path = os.path.join(script_dir, 'bands.dat.gnu')
        band_data = np.loadtxt(data_path)
    except Exception as e:
        print(f"Error reading file 'bands.dat.gnu': {e}")
        sys.exit(1)  # Exit with a non-zero status to indicate an error

    # Assuming the first column is k-points and the rest are energy bands
    k_points = band_data[:, 0]
    energies = band_data[:, 1:]

    plt.figure(figsize=(8, 6))
    for i in range(energies.shape[1]):
        plt.plot(k_points, energies[:, i], color='blue')
    
    plt.axhline(y=0, color='black', linestyle='--')
    plt.xlabel('Wave Vector')
    plt.ylabel('Energy (eV)')
    plt.title('Band Structure')

    # Define the output path in the same directory as the script
    output_path = os.path.join(script_dir, 'band_structure.png')
    
    # Save the plot to the specified file
    plt.savefig(output_path)
    plt.close()  # Close the plot to free up memory
    print(f"Plot saved to {output_path}")

if __name__ == '__main__':
    plot_bands()
    sys.exit(0)  # Exit the script successfully
