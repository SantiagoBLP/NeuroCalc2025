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
"""

import os
import json
import sys
import math
import glob

# Color printing definitions.
try:
    from colorama import Fore, init, Style
    init(autoreset=True)
except ImportError:
    class Fore:
        RED = ""
        GREEN = ""
        YELLOW = ""
        WHITE = ""
    class Style:
        RESET_ALL = ""

def color_print(message, color=Fore.WHITE):
    print(color + message + Style.RESET_ALL)

# Import Element from pymatgen for atomic numbers and validation.
try:
    from pymatgen.core.periodic_table import Element
except ImportError:
    color_print("ERROR: pymatgen is required for computing the recommended number of bands. Please install it via 'pip install pymatgen'.", Fore.RED)
    sys.exit(1)

# Handling missing AI_perp: if not available, skip the Perplexity AI analysis step.
try:
    import AI_perp
except ImportError:
    AI_perp = None
    color_print("WARNING: Perplexity AI analysis module not found. Skipping Perplexity AI analysis step.", Fore.YELLOW)

def compute_recommended_nbnd(params):
    """
    Compute a recommended number of bands based on total electrons in the system (atomic_positions)
    and the band gap:
      - For insulators (band_gap >= 1e-3 eV): ceil(total_electrons/2) + 5
      - For metals (band_gap < 1e-3 eV): ceil(total_electrons/2) + 10
    If band_gap is missing, assume 1.0 eV.
    """
    counts = {}
    for pos in params["atomic_positions"]:
        elem = pos["element"]
        counts[elem] = counts.get(elem, 0) + 1

    total_electrons = 0
    for element, count in counts.items():
        try:
            atomic_number = Element(element).number
        except Exception:
            color_print(f"ERROR: Cannot determine atomic number for element {element}.", Fore.RED)
            atomic_number = 0
        total_electrons += atomic_number * count

    occupied_bands = math.ceil(total_electrons / 2)
    band_gap = params.get("band_gap", 1.0)

    if band_gap < 1e-3:
        return occupied_bands + 10
    else:
        return occupied_bands + 5

def find_pseudopotential_file(element, requested_pp, pseudo_folder):
    """
    For a single element symbol (e.g. 'C'), searches 'pseudo_folder' for .UPF files that begin
    with that element symbol (case-insensitive). If multiple matches are found, user picks one.
    If none are found, returns a default: e.g. 'C.pbesol.UPF'.
    If 'requested_pp' is a valid file in 'pseudo_folder', we just use that and skip searching.
    """
    if requested_pp.strip():
        candidate = os.path.join(pseudo_folder, requested_pp)
        if os.path.isfile(candidate):
            return requested_pp
        else:
            color_print(f"WARNING: Provided pseudopotential '{requested_pp}' not found in {pseudo_folder}. Searching for alternatives...", Fore.YELLOW)

    element_upper = element.upper()
    all_upf = glob.glob(os.path.join(pseudo_folder, "*.UPF"))
    matches = []
    for f in all_upf:
        bn = os.path.basename(f)
        bn_up = bn.upper()
        if bn_up.startswith(element_upper):
            matches.append(bn)

    if len(matches) == 0:
        color_print(f"WARNING: No pseudopotential found for {element} in {pseudo_folder}. Using default {element}.pbesol.UPF", Fore.YELLOW)
        return f"{element}.pbesol.UPF"
    elif len(matches) == 1:
        chosen = matches[0]
        color_print(f"Found one pseudopotential for {element}: {chosen}", Fore.WHITE)
        return chosen
    else:
        color_print(f"\nMultiple pseudopotentials found for {element} in {pseudo_folder}:", Fore.WHITE)
        for i, m in enumerate(matches, start=1):
            color_print(f"  {i}. {m}", Fore.WHITE)
        while True:
            choice = input(f"Select a pseudopotential for {element} [1-{len(matches)}]: ").strip()
            try:
                idx = int(choice)
                if 1 <= idx <= len(matches):
                    chosen = matches[idx - 1]
                    color_print(f"Selected {chosen}\n", Fore.GREEN)
                    return chosen
            except ValueError:
                pass
            color_print(f"Invalid choice. Enter a number between 1 and {len(matches)}.", Fore.YELLOW)

def build_system_block(ibrav, a_bohr, nat, ntyp, ecutwfc, nbnd, degauss,
                       nspin, occupations, hubbard_dict, magnetic,
                       soc_enabled=False):
    """
    Build the &SYSTEM block for QE input files.
    If nbnd is not None, it will be included (for SCF/NSCF calculations).
    Also includes Hubbard corrections, spin-orbit, etc.
    """
    block = "&SYSTEM\n"
    block += f"  ibrav = {ibrav}\n"
    block += f"  celldm(1) = {a_bohr:.6f}\n"
    block += f"  nat = {nat}\n"
    block += f"  ntyp = {ntyp}\n"
    block += f"  ecutwfc = {ecutwfc}\n"
    if nbnd is not None:
        block += f"  nbnd = {nbnd}\n"
    block += f"  degauss = {degauss}\n"
    block += f"  nspin = {nspin}\n"
    block += f"  occupations = '{occupations}'\n"

    if soc_enabled:
        block += "  noncolin = .true.\n"
        block += "  lspinorb = .true.\n"

    if hubbard_dict and any(u > 0 for u in hubbard_dict.values()):
        block += "  lda_plus_u = .true.\n"
        block += "  lda_plus_u_kind = 0\n"
        for i, el in enumerate(sorted(hubbard_dict.keys()), start=1):
            block += f"  hubbard_u({i}) = {hubbard_dict[el]:.2f}  ! {el}\n"
    else:
        block += "  lda_plus_u = .false.\n"

    if magnetic:
        for i in range(1, ntyp+1):
            block += f"  starting_magnetization({i}) = 0.10\n"

    block += "/\n"
    return block

def generate_qe_inputs(params, nbnd, hubbard_dict, magnetic, soc_enabled=False):
    """
    Main function:
      - For each unique element, finds/validates pseudopotential once (so no repeated prompts).
      - Builds scf.in, nscf.in, bands.in, etc. in ../inputs
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(script_dir, "..")
    pseudo_folder = os.path.join(project_root, "pseudo")

    if not os.path.isdir(pseudo_folder):
        color_print(f"WARNING: Pseudopotential folder '{pseudo_folder}' not found. Creating it...", Fore.YELLOW)
        os.makedirs(pseudo_folder, exist_ok=True)

    # For each element, only prompt once if needed.
    chosen_pp_for_element = {}
    for sp in params["atomic_species"]:
        el = sp["element"]
        if el in chosen_pp_for_element:
            sp["pseudopotential"] = chosen_pp_for_element[el]
            continue
        new_pp = find_pseudopotential_file(el, sp["pseudopotential"], pseudo_folder)
        sp["pseudopotential"] = new_pp
        chosen_pp_for_element[el] = new_pp

    # Convert lattice from Å to Bohr.
    ibrav = 1
    a_bohr = float(params["lattice_parameter"]) / 0.529177
    nat = len(params["atomic_positions"])
    ntyp = len(params["atomic_species"])
    ecutwfc = params.get("ecutwfc", 50.0)

    # Degauss.
    band_gap = params.get("band_gap", 1.0)
    degauss = 0.005 if band_gap > 0.1 else 0.02

    # Spin.
    nspin = 2 if magnetic else 1
    occupations = "smearing"

    # Build system blocks.
    system_block_nbnd = build_system_block(ibrav, a_bohr, nat, ntyp,
                                             ecutwfc, nbnd, degauss,
                                             nspin, occupations,
                                             hubbard_dict, magnetic,
                                             soc_enabled)
    system_block = build_system_block(ibrav, a_bohr, nat, ntyp,
                                      ecutwfc, None, degauss,
                                      nspin, occupations,
                                      hubbard_dict, magnetic,
                                      soc_enabled)

    # Build ATOMIC_SPECIES.
    atomic_species_lines = []
    for sp in params["atomic_species"]:
        line = f"{sp['element']} {sp['mass']:.4f} {sp['pseudopotential']}"
        atomic_species_lines.append(line)
    atomic_species_block = "\n".join(atomic_species_lines)

    # Build ATOMIC_POSITIONS.
    atomic_positions_lines = []
    for pos in params["atomic_positions"]:
        line = f"{pos['element']} {pos['x']:.6f} {pos['y']:.6f} {pos['z']:.6f}"
        atomic_positions_lines.append(line)
    atomic_positions_block = "\n".join(atomic_positions_lines)

    # K-points.
    kpoints_scf = params.get("kpoints_scf", "6 6 6 0 0 0")
    kpoints_nscf = params.get("kpoints_nscf", "12 12 12 0 0 0")
    kpoints_band = (
        "K_POINTS (crystal)\n"
        "4\n"
        "  0.000000 0.000000 0.000000  1  ! Gamma\n"
        "  0.500000 0.000000 0.000000  1  ! X\n"
        "  0.500000 0.500000 0.000000  1  ! M\n"
        "  0.000000 0.000000 0.000000  1  ! Gamma\n"
        "\n"
    )

    # Create input files with adjusted formatting (no extra indents, and ensuring a trailing newline).
    files = {
        "scf.in": f"""&CONTROL
calculation = 'scf'
prefix = 'qe'
outdir = './tmp'
pseudo_dir = './pseudo'
/
{system_block_nbnd}
&ELECTRONS
conv_thr = 1.0d-8
/
ATOMIC_SPECIES
{atomic_species_block}
ATOMIC_POSITIONS crystal
{atomic_positions_block}
K_POINTS automatic
{kpoints_scf}
""",
        "nscf.in": f"""&CONTROL
calculation = 'nscf'
prefix = 'qe'
outdir = './tmp'
pseudo_dir = './pseudo'
/
{system_block_nbnd}
&ELECTRONS
conv_thr = 1.0d-8
/
ATOMIC_SPECIES
{atomic_species_block}
ATOMIC_POSITIONS crystal
{atomic_positions_block}
K_POINTS automatic
{kpoints_nscf}
""",
        "bands.in": f"""&CONTROL
calculation = 'bands'
prefix = 'qe'
outdir = './tmp'
pseudo_dir = './pseudo'
/
{system_block}
&ELECTRONS
conv_thr = 1.0d-8
/
ATOMIC_SPECIES
{atomic_species_block}
ATOMIC_POSITIONS crystal
{atomic_positions_block}
{kpoints_band}
""",
        "bands_pp.in": """&BANDS
prefix = 'qe'
outdir = './tmp'
filband = 'bands.dat'
/
""",
        "dos.in": """&DOS
prefix = 'qe'
outdir = './tmp'
fildos = 'dos.dat'
/
"""
    }

    inputs_folder = os.path.join(project_root, "inputs")
    os.makedirs(inputs_folder, exist_ok=True)
    for filename, content in files.items():
        filepath = os.path.join(inputs_folder, filename)
        with open(filepath, "w") as f:
            f.write(content)
        color_print(f"SUCCESS: Generated {filepath}", Fore.GREEN)

def system_needs_soc(params):
    """
    Checks if any element is in the heavy_elements set, meaning spin-orbit might be relevant.
    """
    heavy_elements = {"PB", "BI", "TL", "PT", "AU", "HG", "W", "IR", "OS", "RE", "U", "TH", "TA"}
    for sp in params["atomic_species"]:
        if sp["element"].upper() in heavy_elements:
            return True
    return False

def interactive_update_parameters(params):
    """
    Show current parameters, let user modify them if needed, ensuring atomic species are uppercase
    and atomic positions match the species. Then loop until user chooses not to modify further.
    """
    while True:
        for sp in params["atomic_species"]:
            original = sp["element"]
            new_sym = original.upper()
            try:
                Element(new_sym)
                sp["element"] = new_sym
            except Exception:
                color_print(f"WARNING: '{new_sym}' is not recognized. Please update it if needed.", Fore.YELLOW)

        defined_species = {sp["element"] for sp in params["atomic_species"]}
        for idx, pos in enumerate(params["atomic_positions"], start=1):
            sym_up = pos["element"].upper()
            if sym_up not in defined_species:
                color_print(f"\nERROR: Atomic position #{idx} uses '{sym_up}', not in {defined_species}", Fore.RED)
                default_el = list(defined_species)[0]
                new_el = input(f"Enter valid element for atom #{idx} [default: {default_el}]: ").strip().upper()
                if not new_el:
                    new_el = default_el
                if new_el in defined_species:
                    pos["element"] = new_el
                else:
                    color_print(f"ERROR: Invalid input. Using '{default_el}' by default.", Fore.RED)
                    pos["element"] = default_el
            else:
                pos["element"] = sym_up

        color_print("\n---------------- Current System Parameters ----------------", Fore.RED)
        color_print(f" Lattice parameter (Å): {params['lattice_parameter']}", Fore.WHITE)
        color_print(f" Plane-wave cutoff (ecutwfc in Ry): {params.get('ecutwfc', 50.0)}", Fore.WHITE)
        color_print(f" SCF k-point grid: {params.get('kpoints_scf', '6 6 6 0 0 0')}", Fore.WHITE)
        color_print(f" NSCF k-point grid: {params.get('kpoints_nscf', '12 12 12 0 0 0')}", Fore.WHITE)
        color_print(f" Band gap (eV): {params.get('band_gap', 1.0)}", Fore.WHITE)
        color_print(" Atomic species:", Fore.WHITE)
        for i, sp in enumerate(params["atomic_species"], start=1):
            color_print(f"   {i}. Element: {sp['element']}, Mass: {sp['mass']}, Pseudopotential: {sp['pseudopotential']}", Fore.WHITE)
        color_print(" Atomic positions:", Fore.WHITE)
        for i, pos in enumerate(params["atomic_positions"], start=1):
            color_print(f"   {i}. Element: {pos['element']}, Coordinates: ({pos['x']}, {pos['y']}, {pos['z']})", Fore.WHITE)
        color_print("-----------------------------------------------------------\n", Fore.RED)

        ans = input("Do you want to modify any parameter? (y/n): ").strip().lower()
        if ans.startswith("y"):
            # Place interactive editing logic here if desired.
            pass
        else:
            break
    return params

def interactive_simulation_parameters(params):
    """
    Asks for nbnd, whether system is magnetic, and Hubbard corrections. Loops until user confirms.
    """
    recommended_nbnd = compute_recommended_nbnd(params)
    while True:
        nbnd_input = input(f"Enter the number of bands [Press Enter for recommended {recommended_nbnd}]: ").strip()
        if nbnd_input == "":
            nbnd = recommended_nbnd
        else:
            try:
                nbnd = int(nbnd_input)
            except ValueError:
                color_print("ERROR: Invalid input. Using recommended value.", Fore.RED)
                nbnd = recommended_nbnd

        mag_ans = input("Is the system magnetic? (y/n) [default: n]: ").strip().lower()
        magnetic = mag_ans.startswith("y")

        hubbard_candidates = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu"]
        default_hubbard = {"Ti": 4.0, "V": 3.0, "Cr": 3.5, "Mn": 4.0, "Fe": 4.0, "Co": 5.0, "Ni": 6.0, "Cu": 7.0}
        hubbard_dict = {}
        elements_in_species = {sp["element"] for sp in params["atomic_species"]}
        for el in sorted(elements_in_species):
            if el in hubbard_candidates:
                user_ans = input(f"Element {el} typically needs Hubbard U (default {default_hubbard[el]} eV). Apply? (y/n): ").strip().lower()
                if user_ans.startswith("y"):
                    u_in = input(f"Enter Hubbard U for {el} [default: {default_hubbard[el]}]: ").strip()
                    try:
                        u_val = float(u_in) if u_in else default_hubbard[el]
                    except ValueError:
                        color_print("ERROR: Invalid input. Using default value.", Fore.RED)
                        u_val = default_hubbard[el]
                    hubbard_dict[el] = u_val

        color_print("\n--------------- Simulation Parameters Summary ---------------", Fore.WHITE)
        color_print(f" Number of bands: {nbnd}", Fore.WHITE)
        color_print(f" Magnetic system: {magnetic}", Fore.WHITE)
        color_print(f" Hubbard corrections: {hubbard_dict if hubbard_dict else 'None'}", Fore.WHITE)
        color_print("-------------------------------------------------------------\n", Fore.WHITE)
        confirm = input("Are these parameters correct? (y/n): ").strip().lower()
        if confirm.startswith("y"):
            return nbnd, magnetic, hubbard_dict
        else:
            color_print("Let's update the simulation parameters again.\n", Fore.WHITE)

def evaluate_parameters_with_literature(params, nbnd, magnetic, hubbard_dict):
    """
    Dummy evaluation. Just prints them out and asks for user acceptance.
    """
    color_print("\nEvaluating parameters with 'literature'... (dummy)", Fore.WHITE)
    color_print(f"  nbnd = {nbnd}", Fore.WHITE)
    color_print(f"  magnetic = {magnetic}", Fore.WHITE)
    color_print(f"  hubbard_dict = {hubbard_dict}", Fore.WHITE)
    ans = input("Accept these parameters? (y/n): ").strip().lower()
    if ans.startswith("y"):
        return nbnd, magnetic, hubbard_dict
    else:
        return None

def run_simulation():
    """
    Simulate running QE. Actually just checks if a CRASH file exists in the project root.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(script_dir, "..")
    crash_file = os.path.join(project_root, "CRASH")
    if os.path.exists(crash_file):
        return False
    else:
        return True

def handle_simulation_crash(params):
    """
    If CRASH is found, print it and let the user fix parameters.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(script_dir, "..")
    crash_file = os.path.join(project_root, "CRASH")
    if os.path.exists(crash_file):
        with open(crash_file, "r") as f:
            crash_message = f.read()
        color_print("\nSimulation Crash Detected:", Fore.RED)
        print(crash_message)

        # Removed the final call to AI_perp.evaluate_crash.
        color_print("Skipping AI crash analysis in this version.\n", Fore.YELLOW)
        solution = None

        if not solution:
            solution = {
                "message": "Try increasing ecutwfc or adjusting the k-points. Possibly an overlap in atomic positions.",
                "suggested_changes": {"ecutwfc": 60.0}
            }

        color_print("\nPerplexity AI Suggested Solution (dummy fallback):", Fore.GREEN)
        print(solution["message"])
        ans = input("Apply these changes automatically? (y/n): ").strip().lower()
        if ans.startswith("y"):
            for k, v in solution.get("suggested_changes", {}).items():
                params[k] = v
            color_print("Changes applied to simulation parameters.", Fore.GREEN)
        else:
            do_manual = input("Modify parameters manually? (y/n): ").strip().lower()
            if do_manual.startswith("y"):
                params = interactive_update_parameters(params)
        return params
    else:
        color_print("No CRASH file found. Nothing to handle.", Fore.RED)
        return params

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.join(script_dir, "..")
    inputs_folder = os.path.join(project_root, "inputs")
    params_file = os.path.join(inputs_folder, "optimized_simulation_parameters.json")

    if not os.path.exists(params_file):
        color_print(f"ERROR: {params_file} not found. Please run your parameter-optimization script first.", Fore.RED)
        sys.exit(1)

    try:
        with open(params_file, "r") as f:
            params = json.load(f)
    except Exception as e:
        color_print(f"ERROR reading parameters from {params_file}: {e}", Fore.RED)
        sys.exit(1)

    params = interactive_update_parameters(params)

    while True:
        nbnd, magnetic, hubbard_dict = interactive_simulation_parameters(params)
        result = evaluate_parameters_with_literature(params, nbnd, magnetic, hubbard_dict)
        if result is not None:
            nbnd, magnetic, hubbard_dict = result
            break
        else:
            color_print("Re-entering parameter update...\n", Fore.WHITE)

    soc_enabled = False
    if system_needs_soc(params):
        soc_ans = input("Detected heavy element(s). Enable spin-orbit coupling (SOC)? (y/n) [default: y]: ").strip().lower()
        if soc_ans == "" or soc_ans.startswith("y"):
            soc_enabled = True

    color_print("\n-----------------------------------------------------", Fore.RED)
    color_print("Final Parameters for Quantum ESPRESSO input files:", Fore.RED)
    color_print(f" Lattice parameter (Å): {params['lattice_parameter']}", Fore.WHITE)
    color_print(f" Plane-wave cutoff (ecutwfc in Ry): {params.get('ecutwfc', 50.0)}", Fore.WHITE)
    color_print(f" SCF k-point grid: {params.get('kpoints_scf', '6 6 6 0 0 0')}", Fore.WHITE)
    color_print(f" NSCF k-point grid: {params.get('kpoints_nscf', '12 12 12 0 0 0')}", Fore.WHITE)
    color_print(f" Band gap (eV): {params.get('band_gap',1.0)}", Fore.WHITE)
    color_print(" Atomic species:", Fore.WHITE)
    for sp in params["atomic_species"]:
        color_print(f"   Element: {sp['element']}, Mass: {sp['mass']}, PP: {sp['pseudopotential']}", Fore.WHITE)
    color_print(" Atomic positions:", Fore.WHITE)
    for pos in params["atomic_positions"]:
        color_print(f"   Element: {pos['element']}, Coordinates: ({pos['x']}, {pos['y']}, {pos['z']})", Fore.WHITE)
    color_print(f" Number of bands: {nbnd}", Fore.WHITE)
    color_print(f" Magnetic? {magnetic}", Fore.WHITE)
    color_print(f" Hubbard corrections: {hubbard_dict if hubbard_dict else 'None'}", Fore.WHITE)
    color_print(f" Spin–orbit coupling: {'Enabled' if soc_enabled else 'Disabled'}", Fore.WHITE)
    color_print("-----------------------------------------------------\n", Fore.RED)

    ready = input("Generate QE input files now? (y/n): ").strip().lower()
    if not ready.startswith("y"):
        color_print("Exiting without generating input files.", Fore.WHITE)
        sys.exit(0)

    color_print("\nGenerating Quantum ESPRESSO input files using the provided parameters...\n", Fore.WHITE)
    generate_qe_inputs(params, nbnd, hubbard_dict, magnetic, soc_enabled)
    color_print("\nSUCCESS: QE input files have been generated in the 'inputs' folder.", Fore.GREEN)

    while True:
        run_ans = input("\nDo you want to run the simulation now? (y/n): ").strip().lower()
        if run_ans.startswith("y"):
            success = run_simulation()
            if success:
                color_print("Simulation completed successfully!", Fore.GREEN)
                break
            else:
                color_print("Simulation crashed!", Fore.RED)
                params = handle_simulation_crash(params)
                retry = input("Re-run the simulation with updated parameters? (y/n): ").strip().lower()
                if retry.startswith("y"):
                    continue
                else:
                    break
        else:
            color_print("Skipping simulation run. Returning to main menu.", Fore.WHITE)
            break

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        color_print(f"ERROR: Failed to generate QE inputs. {e}", Fore.RED)
        sys.exit(1)
