#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import sys  # Import the sys module to use sys.exit()

# Constants
UNITS_DATA = 1  # 1 => data in eV already, 0 => Rydbergs => convert to eV
RY_TO_EV = 13.605698

# Edit these for your runs
FILES = {
    'lowq': 'qe.plot_eps.dat',     # Default name after si.spectrum_lowq.in
}

def load_eps_file(filename):
    """
    Loads the typical output from turbo_spectrum.x:
      col0 = Energy, col1 = Re(eps), col2 = Im(eps),
      col3 = Re(1/eps), col4 = Im(1/eps).
    Returns energy(eV), re_eps, im_eps, re_eps_inv, im_eps_inv
    """
    data = np.loadtxt(filename)
    energy_in = data[:, 0]
    re_epsilon = data[:, 1]
    im_epsilon = data[:, 2]
    re_epsilon_inv = data[:, 3]
    im_epsilon_inv = data[:, 4]

    if UNITS_DATA == 0:
        energy = energy_in * RY_TO_EV
    else:
        energy = energy_in

    return energy, re_epsilon, im_epsilon, re_epsilon_inv, im_epsilon_inv

def plot_eels(energy, re_eps, im_eps, re_eps_inv, im_eps_inv, q_label='Q'):
    """
    Plots real & imag parts of eps, 1/eps, and the loss function.
    """
    loss = -im_eps_inv

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()

    # Real eps
    axes[0].plot(energy, re_eps, 'b-')
    axes[0].set_xlabel('Energy (eV)')
    axes[0].set_ylabel('Re[ε(ω)]')
    axes[0].set_title(f'Re[ε(ω)] @ {q_label}')
    axes[0].grid(True)

    # Imag eps
    axes[1].plot(energy, im_eps, 'r-')
    axes[1].set_xlabel('Energy (eV)')
    axes[1].set_ylabel('Im[ε(ω)]')
    axes[1].set_title(f'Im[ε(ω)] @ {q_label}')
    axes[1].grid(True)

    # Real 1/eps
    axes[2].plot(energy, re_eps_inv, 'g-')
    axes[2].set_xlabel('Energy (eV)')
    axes[2].set_ylabel('Re[1/ε(ω)]')
    axes[2].set_title(f'Re[1/ε(ω)] @ {q_label}')
    axes[2].grid(True)

    # Loss = -Im[1/eps]
    axes[3].plot(energy, loss, 'k-')
    axes[3].set_xlabel('Energy (eV)')
    axes[3].set_ylabel('Loss')
    axes[3].set_title(f'EELS Loss = -Im[1/ε(ω)] @ {q_label}')
    axes[3].grid(True)

    fig.tight_layout()
    # Set x-axis limits to 0-100 eV for all subplots
    for ax in axes:
        ax.set_xlim(0, 100)

    # Save the plot to the same directory as the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_filename = f'eels_plot_{q_label}.png'
    output_path = os.path.join(script_dir, output_filename)
    plt.savefig(output_path)
    plt.close()
    print(f"Plot saved to {output_path}")

if __name__ == '__main__':
    # Low-Q data
    try:
        (en_low,
         re_eps_low, im_eps_low,
         re_epsinv_low, im_epsinv_low) = load_eps_file(FILES['lowq'])
        plot_eels(en_low, re_eps_low, im_eps_low,
                  re_epsinv_low, im_epsinv_low, q_label='lowQ')
    except Exception as e:
        print("Could not plot low-q data:", e)

    sys.exit(0)  # Exit the script gracefully
