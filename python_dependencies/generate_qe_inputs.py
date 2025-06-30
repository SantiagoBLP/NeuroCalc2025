#!/usr/bin/env python
import os
import glob
import sys
import math
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# --- Optional Color Support ---
try:
    from colorama import init, Fore, Style
    init(autoreset=True)
except ImportError:
    # If colorama is not available, define dummies so the script still runs
    class Fore:
        RED = ""
        GREEN = ""
        YELLOW = ""
        WHITE = ""
        CYAN = ""
    class Style:
        BRIGHT = ""
        RESET_ALL = ""

def color_print(text, color=Fore.WHITE):
    print(color + text + Style.RESET_ALL)

def center_text(text, width=None):
    """
    Utility to center a string for nicer console output.
    """
    if width is None:
        try:
            width = shutil.get_terminal_size().columns
        except Exception:
            width = 80
    return text.center(width)

def print_header(header):
    """
    Print a header in RED with a centered title.
    """
    line = "=" * 80
    color_print(center_text(line), Fore.RED)
    color_print(center_text(header), Fore.RED)
    color_print(center_text(line), Fore.RED)

def print_divider():
    """
    Print a horizontal divider in plain white.
    """
    color_print(center_text("-" * 80), Fore.WHITE)

#########################################
# Attempt to Import AI_perp for Perplexity AI
#########################################
try:
    import AI_perp
except ImportError:
    AI_perp = None
    color_print("[WARNING] The AI_perp module was not found. The Perplexity AI step will be skipped.", Fore.YELLOW)

#########################################
# Attempt to Import Pymatgen Modules
#########################################
try:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.periodic_table import Element
    # The next lines are used by get_kpath_from_materials_project:
    from mp_api.client import MPRester
    from pymatgen.symmetry.bandstructure import HighSymmKpath
except ImportError:
    color_print("pymatgen is required. Please install it via 'pip install pymatgen'.", Fore.RED)
    sys.exit(1)

#########################################
# PROJECT SETUP: Define paths relative to the project root
#########################################
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# Absolute paths (for file operations)
API_KEY_FILE   = os.path.join(PROJECT_ROOT, "api_keys", "api_key_mp.txt")
PSEUDO_DIR     = os.path.join(PROJECT_ROOT, "pseudo")
TMPDIR         = os.path.join(PROJECT_ROOT, "tmp")
INPUTS_FOLDER  = os.path.join(PROJECT_ROOT, "inputs")

# Ensure that the tmp and inputs folders exist.
os.makedirs(TMPDIR, exist_ok=True)
os.makedirs(INPUTS_FOLDER, exist_ok=True)

if not os.path.isfile(API_KEY_FILE):
    color_print(f"Error: {API_KEY_FILE} not found. Please create the file with your Materials Project API key.", Fore.RED)
    sys.exit(1)
with open(API_KEY_FILE, "r", encoding="utf-8") as f:
    API_KEY = f.read().strip()

#########################################
# Define Relative Paths for QE Input Files
#########################################
TMPDIR_REL     = "./tmp"
PSEUDO_DIR_REL = "./pseudo/"

#########################################
# Helper Function for Pseudopotentials
#########################################
def get_pseudopotential_filename(symbol, for_soc=False):
    """
    Searches PSEUDO_DIR for UPF files matching the given element symbol.
    If for_soc=True, look for a "rel" or "spin-orbit" type of file.
    If none found, fallback to a default. For SOC, we warn if not found.
    """
    if for_soc:
        pattern = os.path.join(PSEUDO_DIR, f"{symbol}*_rel*.UPF")
    else:
        pattern = os.path.join(PSEUDO_DIR, f"{symbol}*.UPF")
    files = glob.glob(pattern)
    
    if not files:
        if for_soc:
            color_print(f"[WARNING] No SOC-appropriate UPF file found for {symbol} in {PSEUDO_DIR}.", Fore.YELLOW)
            default_choice = input(Fore.YELLOW + f"Continue with the non-relativistic pseudopotential for {symbol}? (y/n): ").strip().lower()
            if not default_choice.startswith("y"):
                color_print("Aborting due to missing SOC pseudopotential.", Fore.RED)
                sys.exit(1)
        pseudo_filename = f"{symbol}.pbesol-spn-kjpaw_psl.1.0.0.UPF"
        if symbol == "O":
            pseudo_filename = f"{symbol}.pbesol-n-kjpaw_psl.1.0.0.UPF"
        color_print(f"[INFO] No UPF file found for {symbol} in {PSEUDO_DIR}. Using default: {pseudo_filename}", Fore.WHITE)
        return pseudo_filename

    if len(files) == 1:
        chosen = os.path.basename(files[0])
        color_print(f"[INFO] One UPF file found for {symbol}: {chosen}", Fore.WHITE)
        return chosen

    color_print(f"\nMultiple pseudopotential files found for {symbol}:", Fore.WHITE)
    for idx, filepath in enumerate(files, start=1):
        color_print(f"  {idx}. {os.path.basename(filepath)}", Fore.WHITE)
    while True:
        try:
            choice = int(input("Select the number of the desired pseudopotential: "))
            if 1 <= choice <= len(files):
                chosen = os.path.basename(files[choice - 1])
                color_print(f"Selected: {chosen}\n", Fore.GREEN)
                return chosen
            else:
                color_print("Invalid selection. Please enter a number from the list above.", Fore.YELLOW)
        except ValueError:
            color_print("Invalid input. Please enter a valid number.", Fore.YELLOW)

#########################################
# Retrieve High-Symmetry k–Path from Materials Project
#########################################
def get_kpath_from_materials_project(material_id):
    color_print("🔍 [INFO] Fetching high-symmetry k-path for material {} ...".format(material_id), Fore.WHITE)
    with MPRester(API_KEY) as mpr:
        structure = mpr.get_structure_by_material_id(material_id)
    analyzer = SpacegroupAnalyzer(structure)
    primitive_structure = analyzer.get_primitive_standard_structure()
    kpath_obj = HighSymmKpath(primitive_structure)
    kpoints = kpath_obj.kpath["kpoints"]
    path = kpath_obj.kpath["path"]

    color_print("[INFO] High-symmetry k-points:", Fore.WHITE)
    for label, coords in kpoints.items():
        color_print(f"  {label}: {coords}", Fore.WHITE)
    color_print("[INFO] Path segments:", Fore.WHITE)
    for segment in path:
        color_print("  " + " -> ".join(segment), Fore.WHITE)
    
    # Save k-path reference file
    kpath_file = os.path.join(INPUTS_FOLDER, "kpath_qe.txt")
    with open(kpath_file, "w") as f:
        f.write("K_POINTS (crystal)\n")
        f.write(f"{len(kpoints)}\n")
        for label, coords in kpoints.items():
            f.write(f"  {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}  1  ! {label}\n")
    color_print(f"[INFO] K-path saved to {kpath_file}", Fore.WHITE)
    return kpoints, path

#########################################
# Retrieve Material Data (including extra properties)
#########################################
def fetch_material_data(material_id):
    color_print("[I'm searching trough my database..!] Fetching material data for {} ...".format(material_id), Fore.WHITE)
    with MPRester(API_KEY) as mpr:
        material = mpr.materials.summary.search(
            material_ids=[material_id],
            fields=["material_id", "formula_pretty", "band_gap", "structure", "density", "symmetry", "is_magnetic"]
        )
        if material:
            color_print("[INFO] Material data fetched successfully.", Fore.WHITE)
            return material[0]
        else:
            color_print("No data found for the given Material ID.", Fore.RED)
            sys.exit(1)

#########################################
# Helper Function to Compute Recommended Number of Bands
#########################################
def compute_recommended_nbnd(structure, band_gap):
    total_electrons = 0
    for el, amt in structure.composition.items():
        total_electrons += el.number * amt
    occupied_bands = math.ceil(total_electrons / 2)
    if band_gap < 1e-3:
        recommended = occupied_bands + 10
    else:
        recommended = occupied_bands + 5
    return recommended

#########################################
# Determine ibrav based on crystal symmetry
#########################################
def determine_ibrav(symmetry):
    bravais_lattice_mapping = {
        "cubic": 1,
        "fcc": 2,
        "bcc": 3,
        "hexagonal": 4,
        "trigonal": 5,
        "tetragonal": 6,
        "bct": 7,
        "orthorhombic": 0,
        "base-centered orthorhombic": 0,
        "face-centered orthorhombic": 0,
        "body-centered orthorhombic": 0,
        "monoclinic": 0,
        "base-centered monoclinic": 0,
        "triclinic": 0
    }
    crystal_sys = str(symmetry.crystal_system).lower()
    return bravais_lattice_mapping.get(crystal_sys, 0)

#########################################
# &SYSTEM Block Builder
#########################################
def build_system_block_improved(ibrav, nat, ntyp, structure, ecutwfc, nbnd, degauss,
                                nspin, occupations, hubbard_dict=None, magnetic=False,
                                spin_orbit=False):
    block = "&SYSTEM\n"
    block += f"  ibrav = {ibrav}\n"
    if ibrav != 0:
        if ibrav in [1, 2, 3]:
            a_bohr = structure.lattice.a / 0.529177
            block += f"  celldm(1) = {a_bohr:.6f}\n"
        elif ibrav == 4:
            a_bohr = structure.lattice.a / 0.529177
            c_over_a = structure.lattice.c / structure.lattice.a
            block += f"  celldm(1) = {a_bohr:.6f}\n"
            block += f"  celldm(3) = {c_over_a:.6f}\n"
        elif ibrav == 5:
            a_bohr = structure.lattice.a / 0.529177
            alpha_radians = math.radians(structure.lattice.alpha)
            cos_alpha = math.cos(alpha_radians)
            block += f"  celldm(1) = {a_bohr:.6f}\n"
            block += f"  celldm(2) = {cos_alpha:.6f}\n"
        elif ibrav in [6, 7]:
            a_bohr = structure.lattice.a / 0.529177
            c_over_a = structure.lattice.c / structure.lattice.a
            block += f"  celldm(1) = {a_bohr:.6f}\n"
            block += f"  celldm(3) = {c_over_a:.6f}\n"
        else:
            a_bohr = structure.lattice.a / 0.529177
            block += f"  celldm(1) = {a_bohr:.6f}\n"
    block += f"  nat = {nat}\n"
    block += f"  ntyp = {ntyp}\n"
    block += f"  ecutwfc = {ecutwfc}\n"
    block += f"  nbnd = {nbnd}\n"
    block += f"  degauss = {degauss}\n"
    block += f"  nspin = {nspin}\n"
    block += f"  occupations = '{occupations}'\n"
    if spin_orbit:
        block += f"  noncolin = .true.\n"
        block += f"  lspinorb = .true.\n"
    if hubbard_dict and any(val > 0 for val in hubbard_dict.values()):
        block += f"  lda_plus_u = .true.\n"
        block += f"  lda_plus_u_kind = 0\n"
        for i, el in enumerate(sorted(hubbard_dict.keys()), start=1):
            u_val = hubbard_dict[el]
            block += f"  hubbard_u({i}) = {u_val:.2f}  ! {el}\n"
    else:
        block += f"  lda_plus_u = .false.\n"
    if magnetic:
        for i in range(1, ntyp + 1):
            block += f"  starting_magnetization({i}) = 0.10\n"
    block += "/\n"
    return block

#########################################
# Helper to validate k-point grid input
#########################################
def validate_kgrid(kgrid_input, default):
    parts = kgrid_input.split()
    if len(parts) == 3:
        return kgrid_input + " 0 0 0"
    elif len(parts) == 6:
        return kgrid_input
    else:
        color_print("Invalid k-point grid input. Using default.", Fore.YELLOW)
        default_parts = default.split()
        if len(default_parts) == 3:
            return default + " 0 0 0"
        else:
            return default

#########################################
# QE Input File Generation
#########################################
def generate_qe_inputs(material_data, nbnd, kgrid_scf, kgrid_nscf, hubbard_dict, spin_orbit):
    color_print("[I'm thinking...] Generating Quantum ESPRESSO input files...", Fore.WHITE)
    ibrav = determine_ibrav(material_data.symmetry)
    
    # Heuristic for ecutwfc and degauss
    if material_data.band_gap > 0.1:
        ecutwfc = 80.0
        degauss = 0.005
    else:
        ecutwfc = 50.0
        degauss = 0.02

    # Magnetic or not
    if hasattr(material_data, "is_magnetic") and material_data.is_magnetic:
        nspin = 2
        magnetic = True
    else:
        nspin = 1
        magnetic = False

    occupations = "smearing"

    # ATOMIC_SPECIES
    # Use each element's proper symbol (e.g., "Sr", "Ti") as provided by pymatgen.
    elements = sorted({el.symbol for el in material_data.structure.composition.elements})
    atomic_species_lines = []
    for el in elements:
        try:
            mass = Element(el).atomic_mass
        except Exception as e:
            color_print(f"[WARNING] Element {el} not recognized. Please input the correct element symbol:", Fore.YELLOW)
            corrected = input().strip()
            try:
                mass = Element(corrected).atomic_mass
                el = corrected  # Use the corrected symbol as provided
            except Exception as e:
                color_print(f"ERROR: Provided element {corrected} is still not recognized. Aborting.", Fore.RED)
                sys.exit(1)
        pseudo_file = get_pseudopotential_filename(el, for_soc=spin_orbit)
        atomic_species_lines.append(f"{el} {mass:.4f} {pseudo_file}")
    atomic_species = "\n".join(atomic_species_lines)

    # ATOMIC_POSITIONS
    # Use the element's symbol with proper capitalization as provided by the structure.
    atomic_positions = "\n".join(
        f"  {site.specie.symbol} {site.frac_coords[0]:.6f} {site.frac_coords[1]:.6f} {site.frac_coords[2]:.6f}"
        for site in material_data.structure.sites
    )
    
   
    nat = len(material_data.structure.sites)
    ntyp = len(elements)
    if ibrav == 0:
        lattice = material_data.structure.lattice.matrix
        cell_params = "CELL_PARAMETERS angstrom\n" + "\n".join(
            "  " + " ".join(f"{x:.6f}" for x in vec) for vec in lattice
        )
    else:
        cell_params = ""
    
    system_block = build_system_block_improved(ibrav, nat, ntyp, material_data.structure,
                                               ecutwfc, nbnd, degauss, nspin, occupations,
                                               hubbard_dict=hubbard_dict,
                                               magnetic=magnetic,
                                               spin_orbit=spin_orbit)
    
    kpoints_scf  = f"K_POINTS automatic\n  {kgrid_scf}\n"
    kpoints_nscf = f"K_POINTS automatic\n  {kgrid_nscf}\n"
    
    # High-symmetry k-path
    mp_kpoints, mp_path = get_kpath_from_materials_project(material_data.material_id)
    kpoints_band_block = "K_POINTS (crystal)\n"
    kpoints_band_block += f"{len(mp_kpoints)}\n"
    for label, coords in mp_kpoints.items():
        kpoints_band_block += "  " + " ".join(f"{coord:.6f}" for coord in coords) + "  1  ! " + label + "\n"
    
    # Add timer message before generating files
    color_print("⏳ Working in the background, give me one sec...", Fore.GREEN)
    time.sleep(1)
    
    files = {
        "scf.in": f"""&CONTROL
  calculation = 'scf'
  prefix = '{material_data.formula_pretty.replace(' ', '_')}'
  outdir = '{TMPDIR_REL}'
  pseudo_dir = '{PSEUDO_DIR_REL}'
/
{system_block}
&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.7
/
ATOMIC_SPECIES
{atomic_species}
ATOMIC_POSITIONS crystal
{atomic_positions}
{cell_params}
{kpoints_scf}
""",
        "nscf.in": f"""&CONTROL
  calculation = 'nscf'
  prefix = '{material_data.formula_pretty.replace(' ', '_')}'
  outdir = '{TMPDIR_REL}'
  pseudo_dir = '{PSEUDO_DIR_REL}'
/
{system_block}
&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.7
/
ATOMIC_SPECIES
{atomic_species}
ATOMIC_POSITIONS crystal
{atomic_positions}
{cell_params}
{kpoints_nscf}
""",
        "bands.in": f"""&CONTROL
  calculation = 'bands'
  prefix = '{material_data.formula_pretty.replace(' ', '_')}'
  outdir = '{TMPDIR_REL}'
  pseudo_dir = '{PSEUDO_DIR_REL}'
/
{system_block}
&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.7
/
ATOMIC_SPECIES
{atomic_species}
ATOMIC_POSITIONS crystal
{atomic_positions}
{cell_params}
{kpoints_band_block}
""",
        "bands_pp.in": f"""&BANDS
  prefix = '{material_data.formula_pretty.replace(' ', '_')}'
  outdir = '{TMPDIR_REL}'
  filband = 'bands.dat'
/
""",
        "dos.in": f"""&DOS
  prefix = '{material_data.formula_pretty.replace(' ', '_')}'
  outdir = '{TMPDIR_REL}'
  fildos = 'dos.dat'
/
"""
    }
    
    for filename, content in files.items():
        filepath = os.path.join(INPUTS_FOLDER, filename)
        with open(filepath, "w") as f:
            f.write(content)
        color_print(f"[INFO] Wrote {filename} to {filepath}", Fore.WHITE)
    
    color_print("[SUCCESS] All Quantum ESPRESSO input files generated successfully.", Fore.GREEN)
    return os.path.join(INPUTS_FOLDER, "scf.in")


#########################################
# Helper to coloring AI response
#########################################
def print_ai_suggestions_in_colors(ai_text):
    """
    Prints AI suggestions with color-coding:
      - Lines containing 'Suggestion' or starting with '•' in green
      - Lines containing 'DOI', 'ACS Nano', or 'References' in light blue
      - Everything else in white
    """
    lines = ai_text.splitlines()
    for line in lines:
        stripped_line = line.strip().lower()

        # Check if the line is referencing references/DOIs
        if "doi" in stripped_line or "acs nano" in stripped_line or "references" in stripped_line:
            color_print(line, Fore.CYAN)  # Light blue
        # Check if the line is a suggestion
        elif stripped_line.startswith("suggestion") or stripped_line.startswith("•"):
            color_print(line, Fore.GREEN)
        else:
            color_print(line, Fore.WHITE)



#########################################
# Use Perplexity AI to compare generated input with best practices
#########################################
def analyze_input_with_perplexity(input_filepath):
    if AI_perp is None:
        color_print("[WARNING] Skipping AI analysis because AI_perp is missing.", Fore.YELLOW)
        return
    try:
        with open(input_filepath, "r") as f:
            content = f.read()
    except Exception as e:
        color_print(f"Failed to read {input_filepath}: {e}", Fore.RED)
        return

    prompt = (
        "You are an expert in computational materials science. Analyze the following Quantum ESPRESSO SCF input file "
        "and compare it with the best practices from recent literature. Provide detailed suggestions for improvement "
        "with full ACS Nano style references (including DOI). Here is the input file:\n\n"
        f"{content}\n\n"
        "Your analysis:"
    )
    print_divider()
    color_print(center_text("AI ANALYSIS"), Fore.RED)
    print_divider()

    print("I'll check out your code and help you out with the best tips and practices.")
    print("Give me just a sec... Thinking and searching—it won’t take long!")

    # Query the AI
    response = AI_perp.query_perplexity(prompt, max_tokens=4096)

    # Post-process the AI response for color-coded suggestions & references
    print_ai_suggestions_in_colors(response)

    print_divider()
    color_print(center_text("REFERENCES (see above)"), Fore.CYAN)
    print_divider()



#########################################
# Interactive Parameter Tuning Loop
#########################################
def interactive_parameter_tuning(material_data, recommended_nbnd):
    # Prompt explicitly for number of bands first
    user_nb = input(Fore.GREEN + f"Enter the number of bands (default {recommended_nbnd}): ").strip()
    if user_nb:
        try:
            nbnd = int(user_nb)
        except ValueError:
            color_print("Invalid input, using recommended value.", Fore.YELLOW)
            nbnd = recommended_nbnd
    else:
        nbnd = recommended_nbnd

    default_kgrid_scf = "6 6 6"
    kgrid_scf = validate_kgrid(default_kgrid_scf, default_kgrid_scf)
    default_kgrid_nscf = "12 12 12"
    kgrid_nscf = validate_kgrid(default_kgrid_nscf, default_kgrid_nscf)
    default_soc = False
    spin_orbit = default_soc

    # Hubbard references
    hubbard_candidates = ["TI", "V", "CR", "MN", "FE", "CO", "NI", "CU"]
    default_hubbard = {"TI": 4.0, "V": 3.0, "CR": 3.5, "MN": 4.0, "FE": 4.0, "CO": 5.0, "NI": 6.0, "CU": 7.0}
    elements = sorted({str(el).upper() for el in material_data.structure.composition.elements})
    hubbard_dict = {}
    for el in elements:
        if el in hubbard_candidates:
            ans = input(Fore.GREEN + f"✏️  Element {el} typically requires Hubbard corrections (default U = {default_hubbard[el]} eV). Apply Hubbard correction for {el}? (y/n): ").strip().lower()
            if ans.startswith("y"):
                u_input = input(Fore.GREEN + f"✏️  Enter Hubbard U value for {el} [default: {default_hubbard[el]}]: ").strip()
                try:
                    u_val = float(u_input) if u_input else default_hubbard[el]
                except ValueError:
                    color_print("Invalid input. Using default value.", Fore.YELLOW)
                    u_val = default_hubbard[el]
                hubbard_dict[el] = u_val

    while True:
        print_header("CURRENT PARAMETERS")
        color_print(f"1. Number of bands (nbnd): {nbnd}", Fore.GREEN)
        color_print(f"2. SCF k-point grid: {kgrid_scf}", Fore.GREEN)
        color_print(f"3. NSCF k-point grid: {kgrid_nscf}", Fore.GREEN)
        color_print("4. Hubbard corrections:", Fore.GREEN)
        if hubbard_dict:
            for el, u in hubbard_dict.items():
                color_print(f"   {el}: {u} eV", Fore.GREEN)
        else:
            color_print("   None", Fore.GREEN)
        color_print(f"5. Spin Orbit Coupling (SOC): {'Enabled' if spin_orbit else 'Disabled'}", Fore.GREEN)

        ans = input(Fore.GREEN + "Are you satisfied with these parameters? (y/n): ").strip().lower()
        if ans.startswith("y"):
            break
        else:
            param_to_change = input(Fore.GREEN + "Enter the number(s) of the parameter(s) you want to change (e.g., 1 or 2,3,5): ").strip()
            choices = [s.strip() for s in param_to_change.split(",")]
            for choice in choices:
                if choice == "1":
                    nbnd_input = input(Fore.GREEN + f"✏️  Enter the number of bands [default: {recommended_nbnd}]: ").strip()
                    if nbnd_input == "":
                        nbnd = recommended_nbnd
                    else:
                        try:
                            nbnd = int(nbnd_input)
                        except ValueError:
                            color_print("Invalid input. Using default recommended value.", Fore.YELLOW)
                            nbnd = recommended_nbnd
                elif choice == "2":
                    new_kgrid_scf = input(Fore.GREEN + f"✏️  Enter the SCF k-point grid [default: {default_kgrid_scf}]: ").strip()
                    kgrid_scf = validate_kgrid(new_kgrid_scf if new_kgrid_scf else default_kgrid_scf, default_kgrid_scf)
                elif choice == "3":
                    new_kgrid_nscf = input(Fore.GREEN + f"✏️  Enter the NSCF k-point grid [default: {default_kgrid_nscf}]: ").strip()
                    kgrid_nscf = validate_kgrid(new_kgrid_nscf if new_kgrid_nscf else default_kgrid_nscf, default_kgrid_nscf)
                elif choice == "4":
                    new_hubbard_dict = {}
                    for el in elements:
                        if el in hubbard_candidates:
                            ans = input(Fore.GREEN + f"✏️  Element {el} typically requires Hubbard corrections (default U = {default_hubbard[el]} eV). Apply correction? (y/n): ").strip().lower()
                            if ans.startswith("y"):
                                u_input = input(Fore.GREEN + f"✏️  Enter Hubbard U value for {el} [default: {default_hubbard[el]}]: ").strip()
                                try:
                                    u_val = float(u_input) if u_input else default_hubbard[el]
                                except ValueError:
                                    color_print("Invalid input. Using default value.", Fore.YELLOW)
                                    u_val = default_hubbard[el]
                                new_hubbard_dict[el] = u_val
                    hubbard_dict = new_hubbard_dict
                elif choice == "5":
                    soc_ans = input(Fore.GREEN + "✏️  Enable spin orbit coupling (SOC)? (y/n): ").strip().lower()
                    spin_orbit = soc_ans.startswith("y")
                else:
                    color_print(f"Invalid choice: {choice}", Fore.YELLOW)
    return nbnd, kgrid_scf, kgrid_nscf, hubbard_dict, spin_orbit

#########################################
# Main Routine
#########################################
def main():
    print_header("Starting generate_qe_inputs.py")
    if len(sys.argv) < 2:
        default_material = "mp-149"
        color_print(f"[WARNING] No material_id provided. Using default: {default_material}", Fore.YELLOW)
        material_id = default_material
    else:
        material_id = sys.argv[1]
        color_print(f"[INFO] Using provided material_id: {material_id}", Fore.WHITE)
    
    material_data = fetch_material_data(material_id)
    color_print(f"[INFO] Material: {material_data.formula_pretty} (ID: {material_data.material_id})", Fore.WHITE)
    color_print(f"[INFO] Band Gap: {material_data.band_gap} eV, Density: {material_data.density} g/cm³", Fore.WHITE)
    if hasattr(material_data, "is_magnetic"):
        color_print(f"[INFO] is_magnetic: {material_data.is_magnetic}", Fore.WHITE)
    
    recommended_nbnd = compute_recommended_nbnd(material_data.structure, material_data.band_gap)
    color_print(f"[INFO] Recommended number of bands: {recommended_nbnd}", Fore.WHITE)
    
    # Loop: generate inputs, analyze with AI, then decide to change parameters or run simulation
    while True:
        nbnd, kgrid_scf, kgrid_nscf, hubbard_dict, spin_orbit = interactive_parameter_tuning(material_data, recommended_nbnd)
        scf_filepath = generate_qe_inputs(material_data, nbnd, kgrid_scf, kgrid_nscf, hubbard_dict, spin_orbit)
        
        analyze_input_with_perplexity(scf_filepath)
        
        decision = input(Fore.GREEN + "Based on the AI suggestions, type 'change' to modify parameters or 'run' to proceed with simulation: ").strip().lower()
        if decision.startswith("run"):
            break
        elif decision.startswith("change"):
            color_print("Let's update the parameters based on the suggestions.\n", Fore.YELLOW)
        else:
            color_print("Invalid option. Please type 'change' or 'run'.", Fore.YELLOW)
    
    print_header("Parameters finalized. Proceeding with simulation setup...")
    color_print("[SUCCESS] Script completed successfully.", Fore.GREEN)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        color_print(f"ERROR: Failed to generate QE inputs. {e}", Fore.RED)
        sys.exit(1)
