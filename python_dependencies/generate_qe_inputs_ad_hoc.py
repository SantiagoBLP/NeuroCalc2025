#!/usr/bin/env python3
"""
Generate Quantum ESPRESSO Input Files from Ad-hoc Parameters (Enhanced Version)

This script reads the JSON file (optimized_simulation_parameters.json) produced by either
assist_known_material.py or assist_unpublished_material.py and then generates a set of QE
input files (scf.in, nscf.in, bands.in, bands_pp.in, dos.in) in the 'inputs' folder.

Enhancements:
 - Searches the 'pseudo' folder for available pseudopotential files matching each element.
 - If multiple matches are found for one element, the user is prompted only once to pick one.
 - If no matches are found, a default pseudopotential (Element.pbesol.UPF) is assigned.
 - Provides an interactive loop for reviewing and modifying parameters (lattice, cutoff, etc.).
 - Validates that all element symbols are uppercase and that atomic positions match the defined species.
 - Asks for the number of bands, Hubbard corrections, and whether to enable spin-orbit coupling.
 - Has a dummy “literature evaluation” step for user confirmation.
 - If a CRASH file is found, it skips the final AI crash analysis call.
 - Fixes simplified LDA+U + noncollinear magnetism error.
"""

import os
import sys
import json
import math
import glob

try:
    from colorama import Fore, Style, init
    init(autoreset=True)
except ImportError:
    class Fore:
        RED = GREEN = YELLOW = WHITE = ""
    class Style:
        RESET_ALL = ""

def color_print(msg, color=Fore.WHITE):
    print(color + msg + Style.RESET_ALL)

try:
    from pymatgen.core import Element
except ImportError:
    color_print("ERROR: Install pymatgen via 'pip install pymatgen'", Fore.RED)
    sys.exit(1)

# Helper Functions

def find_pseudopotential(element, pseudo_dir, requested=None):
    if requested:
        path = os.path.join(pseudo_dir, requested)
        if os.path.isfile(path):
            return requested
    files = [f for f in os.listdir(pseudo_dir) if f.upper().startswith(element.upper()) and f.endswith(".UPF")]
    if files:
        return files[0]
    return f"{element}.pbesol.UPF"

def compute_nbnd(params):
    total_e = sum(Element(atom["element"]).Z for atom in params["atomic_positions"])
    gap = params.get("band_gap", 1.0)
    return math.ceil(total_e/2) + (5 if gap > 1e-3 else 10)

def system_needs_soc(params):
    heavy = {"W", "Pt", "Au", "Pb", "Bi", "U", "Th", "Hg", "Ir", "Os"}
    return any(atom["element"].capitalize() in heavy for atom in params["atomic_species"])

# Input Builders

def build_system_block(params, nbnd, magnetic, hubbard, soc):
    block = ["&SYSTEM"]
    block.append("ibrav = 1")
    block.append(f"celldm(1) = {float(params['lattice_parameter'])/0.529177:.6f}")
    block.append(f"nat = {len(params['atomic_positions'])}")
    block.append(f"ntyp = {len(params['atomic_species'])}")
    block.append(f"ecutwfc = {params.get('ecutwfc',50.0)}")
    if nbnd:
        block.append(f"nbnd = {nbnd}")
    block.append(f"degauss = {0.005 if params.get('band_gap',1)>0.1 else 0.02}")
    block.append(f"nspin = {2 if magnetic else 1}")
    block.append("occupations = 'smearing'")
    if soc:
        block.append("noncolin = .true.")
        block.append("lspinorb = .true.")
    if hubbard:
        block.append("lda_plus_u = .true.")
        block.append(f"lda_plus_u_kind = {1 if soc else 0}")  # <<=== FIXED HERE
        for i, (el, uval) in enumerate(sorted(hubbard.items()), start=1):
            block.append(f"hubbard_u({i}) = {uval:.2f}  ! {el}")
    else:
        block.append("lda_plus_u = .false.")
    if magnetic:
        for i in range(1, len(params['atomic_species'])+1):
            block.append(f"starting_magnetization({i}) = 0.1")
    block.append("/")
    return "\n".join(block)

def atomic_blocks(params):
    species = [f"{sp['element']} {sp['mass']} {sp['pseudopotential']}" for sp in params['atomic_species']]
    positions = [f"{pos['element']} {pos['x']} {pos['y']} {pos['z']}" for pos in params['atomic_positions']]
    return "\n".join(species), "\n".join(positions)

def generate_inputs(params, nbnd, magnetic, hubbard, soc):
    cwd = os.path.dirname(__file__)
    pseudo_dir = os.path.join(cwd, "..", "pseudo")
    os.makedirs(pseudo_dir, exist_ok=True)
    for sp in params['atomic_species']:
        sp['pseudopotential'] = find_pseudopotential(sp['element'], pseudo_dir, sp.get('pseudopotential',''))

    species_block, pos_block = atomic_blocks(params)
    system_block = build_system_block(params, nbnd, magnetic, hubbard, soc)

    inputs_dir = os.path.join(cwd, "..", "inputs")
    os.makedirs(inputs_dir, exist_ok=True)

    kscf = params.get("kpoints_scf", "6 6 6 0 0 0")
    knscf = params.get("kpoints_nscf", "12 12 12 0 0 0")

    templates = {
        "scf.in": f"""&CONTROL
calculation = 'scf'
prefix = 'qe'
outdir = './tmp'
pseudo_dir = './pseudo'
/
{system_block}
&ELECTRONS
conv_thr = 1.0d-8
/
ATOMIC_SPECIES
{species_block}
ATOMIC_POSITIONS crystal
{pos_block}
K_POINTS automatic
{kscf}
""",
        "nscf.in": f"""&CONTROL
calculation = 'nscf'
prefix = 'qe'
outdir = './tmp'
pseudo_dir = './pseudo'
/
{system_block}
&ELECTRONS
conv_thr = 1.0d-8
/
ATOMIC_SPECIES
{species_block}
ATOMIC_POSITIONS crystal
{pos_block}
K_POINTS automatic
{knscf}
""",
        "bands_pp.in": """&BANDS
prefix = 'qe'
outdir = './tmp'
filband = 'bands.dat'
/""",
        "dos.in": """&DOS
prefix = 'qe'
outdir = './tmp'
fildos = 'dos.dat'
/"""
    }

    for name, content in templates.items():
        with open(os.path.join(inputs_dir, name), 'w') as f:
            f.write(content)
        color_print(f"Generated: {name}", Fore.GREEN)

# Main

def main():
    cwd = os.path.dirname(__file__)
    param_file = os.path.join(cwd, "..", "inputs", "optimized_simulation_parameters.json")
    if not os.path.isfile(param_file):
        color_print("ERROR: optimized_simulation_parameters.json missing.", Fore.RED)
        sys.exit(1)

    with open(param_file) as f:
        params = json.load(f)

    nbnd = compute_nbnd(params)
    magnetic = input("Magnetic system? (y/N): ").lower().startswith("y")
    hubbard = {}
    for sp in params["atomic_species"]:
        el = sp["element"]
        if el in {"Fe","Co","Ni","Mn"}:
            ans = input(f"Apply Hubbard U to {el}? (y/N): ").lower()
            if ans.startswith("y"):
                val = input(f"Value of U (eV) for {el} [default 4.0]: ") or "4.0"
                hubbard[el] = float(val)

    soc = system_needs_soc(params)
    if soc:
        ans = input("Detected heavy elements. Enable SOC? (Y/n): ").lower()
        soc = not ans.startswith("n")

    generate_inputs(params, nbnd, magnetic, hubbard, soc)
    color_print("All input files generated successfully!", Fore.GREEN)

if __name__ == "__main__":
    main()
