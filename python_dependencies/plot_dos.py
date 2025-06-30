import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def plot_dos():
    """Plot density of states from Quantum Espresso output and save to a file."""
    try:
        # Determine the directory where the script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Load DOS data from file 'dos.dat' located in the same directory
        data_path = os.path.join(script_dir, 'dos.dat')
        dos_data = np.loadtxt(data_path)
    except Exception as e:
        print(f"Error reading file 'dos.dat': {e}")
        sys.exit(1)  # Exit with a non-zero status to indicate an error

    # Ensure data has at least two columns (Energy and DOS)
    if dos_data.shape[1] < 2:
        print("Error: 'dos.dat' should have at least two columns (Energy and DOS).")
        sys.exit(1)  # Exit with a non-zero status to indicate an error

    # Extract energy (x-axis) and DOS (y-axis)
    energy = dos_data[:, 0]
    density_of_states = dos_data[:, 1]

    # Filter data to keep only points within -2 eV to 2 eV
    mask = (energy >= -2) & (energy <= 2)
    energy_filtered = energy[mask]
    dos_filtered = density_of_states[mask]

    # Create figure
    plt.figure(figsize=(8, 6))
    plt.plot(energy_filtered, dos_filtered, color='blue')

    # Add a vertical line at the Fermi level (0 eV)
    plt.axvline(x=0, color='black', linestyle='--')

    plt.xlabel('Energy (eV)')
    plt.ylabel('Density of States')
    plt.title('Density of States (Energy range: -2 to 2 eV)')
    plt.xlim(-2, 2)  # Ensure the x-axis is strictly between -2 and 2 eV

    # Define the output path in the same directory as the script
    output_path = os.path.join(script_dir, 'dos_plot.png')

    # Save the plot
    plt.savefig(output_path)
    plt.close()  # Close the plot to free up memory
    print(f"Plot saved to {output_path}")

if __name__ == '__main__':
    plot_dos()
    sys.exit(0)  # Exit the script successfully
