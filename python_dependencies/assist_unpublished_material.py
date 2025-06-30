#!/usr/bin/env python3
"""
Expert-Level Simulation Parameter Optimization

This script facilitates an expert-level dialogue between a simulation specialist and an AI assistant,
aimed at optimizing parameters for crystal simulations. The assistant engages with you—a distinguished
crystal simulation expert—by systematically collecting key simulation parameters and offering precise,
formal, and technically relevant feedback for optimization. AI recommendations include established research
citations in full ACS Nano style with DOIs and ODI where applicable.

After each evaluation, you may refine the inputs. Once satisfied, the final parameter set is exported as a
structured JSON file ("optimized_simulation_parameters.json") within the inputs folder for subsequent simulation workflows.

Additional behavior in this version:
  • Removes any previously saved JSON before writing a new one.
  • Re-loads the JSON immediately after saving to confirm correctness.
  • Builds the "atomic_species" field as a list of unique species. If a pseudopotential is not provided, a default is assigned.
  • Uses a balanced color scheme: red for headers, yellow for warnings, white for normal text.
"""

import json
import os
import sys
import math
import textwrap
from openai import OpenAI

# Import Element from pymatgen (optional, used for electron count)
try:
    from pymatgen.core.periodic_table import Element
except ImportError:
    print("pymatgen is required. Please install it via 'pip install pymatgen'.")
    sys.exit(1)

# --- Color Setup ---
try:
    from colorama import init, Fore, Style
    init(autoreset=True)
except ImportError:
    class Fore:
        RED = ""
        GREEN = ""
        CYAN = ""
        WHITE = ""
        YELLOW = ""
    class Style:
        BRIGHT = ""
        RESET_ALL = ""

def color_print(text, color=Fore.WHITE):
    print(color + text + Style.RESET_ALL)

def center_text(text, width=None):
    try:
        width = os.get_terminal_size().columns
    except Exception:
        width = 80
    return text.center(width)

def print_header(header):
    line = "=" * 80
    color_print(center_text(line), Fore.RED)
    color_print(center_text(header), Fore.RED)
    color_print(center_text(line), Fore.RED)

def print_divider():
    color_print(center_text("-" * 80), Fore.WHITE)

def print_ai_response_with_references(ai_text):
    """
    Prints the AI response line-by-line:
      - Lines containing "doi", "acs nano", or "references" are printed in cyan;
      - All other lines are printed in white.
    Long lines are wrapped to a generous width.
    """
    wrap_width = 200
    for line in ai_text.splitlines():
        wrapped_line = textwrap.fill(line, width=wrap_width)
        if ("doi" in line.lower()) or ("acs nano" in line.lower()) or ("references" in line.lower()):
            color_print(wrapped_line, Fore.CYAN)
        else:
            color_print(wrapped_line, Fore.WHITE)

def compute_recommended_nbnd(params):
    """
    Compute a recommended number of bands based on the total number of electrons
    in the system (from atomic_positions) and the band gap.
      - For insulators (band_gap >= 1e-3 eV): recommended = ceil(total_electrons/2) + 5.
      - For metals (band_gap < 1e-3 eV): recommended = ceil(total_electrons/2) + 10.
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
            color_print(f"Warning: Cannot determine atomic number for element '{element}'.", Fore.YELLOW)
            atomic_number = 0
        total_electrons += atomic_number * count
    occupied_bands = math.ceil(total_electrons / 2)
    band_gap = params.get("band_gap", 1.0)
    if band_gap < 1e-3:
        return occupied_bands + 10
    else:
        return occupied_bands + 5

def query_perplexity(prompt, max_tokens=512, temperature=0.7, stream=False):
    messages = [
        {"role": "system", "content": "You are an AI assistant specialized in computational materials science. Provide precise, formal, and technically relevant feedback."},
        {"role": "user", "content": prompt}
    ]
    try:
        response = client.chat.completions.create(
            model="sonar-pro",
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
            stream=stream
        )
        if not stream:
            return response.choices[0].message.content.strip()
        else:
            full_response = ""
            for part in response:
                chunk = part.choices[0].delta.content or ""
                full_response += chunk
                print(chunk, end="")
            return full_response
    except Exception as e:
        return f"API request failed: {e}"

# --- Perplexity API Setup ---
API_KEY_FILE = "./api_keys/api_key_perp.txt"
if not os.path.isfile(API_KEY_FILE):
    color_print(f"Error: {API_KEY_FILE} not found. Please create the file with your Perplexity API key.", Fore.RED)
    sys.exit(1)
with open(API_KEY_FILE, "r", encoding="utf-8") as f:
    api_key = f.read().strip()
client = OpenAI(api_key=api_key, base_url="https://api.perplexity.ai")

def print_final_parameters(params):
    print_header("Finalized Simulation Parameters")
    idx = 1
    for key, value in params.items():
        color_print(f"{idx}. {key}: {value}", Fore.WHITE)
        idx += 1
    print_divider()

def main():
    print_header("Expert-Level Simulation Parameter Optimization")
    color_print("Welcome! This tool refines your crystal simulation parameters with AI guidance.\n", Fore.WHITE)

    # Determine project root and inputs folder
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    inputs_folder = os.path.join(project_root, "inputs")
    os.makedirs(inputs_folder, exist_ok=True)
    params_file = os.path.join(inputs_folder, "optimized_simulation_parameters.json")

    # Remove old JSON if it exists
    if os.path.exists(params_file):
        os.remove(params_file)
        color_print(f"Removed old JSON file: {params_file}", Fore.YELLOW)

    # Collect new parameters.
    params = {}
    params["simulation_type"] = "expert"
    params["material_name"] = input("Enter the material name or identifier: ").strip()
    params["material_description"] = input("Enter a description of the material's characteristics (optional): ").strip()
    params["important_properties"] = input("Enter any important properties (optional): ").strip()

    while True:
        try:
            params["lattice_parameter"] = float(input("Enter the lattice parameter (in Å): ").strip())
            params["ecutwfc"] = float(input("Enter the plane-wave cutoff energy (ecutwfc in Ry): ").strip())
            break
        except ValueError:
            color_print("Invalid numerical input. Please re-enter the values.", Fore.YELLOW)

    params["kpoints_scf"] = input("Enter the SCF k-point grid (e.g., '6 6 6 0 0 0'): ").strip()
    params["kpoints_nscf"] = input("Enter the NSCF k-point grid (e.g., '12 12 12 0 0 0'): ").strip()

    try:
        num_species = int(input("Enter the number of distinct atomic species: ").strip())
    except ValueError:
        color_print("Invalid number for atomic species. Exiting.", Fore.RED)
        sys.exit(1)

    atomic_species = []
    for i in range(num_species):
        color_print(f"\nEnter details for atomic species {i+1}:", Fore.WHITE)
        elem = input("  Element symbol (e.g., Fe): ").strip().upper()
        try:
            mass = float(input("  Atomic mass (amu): ").strip())
        except ValueError:
            color_print("Invalid mass input. Exiting.", Fore.RED)
            sys.exit(1)
        pseudo = input("  Pseudopotential filename (e.g., Fe.pbesol.UPF): ").strip()
        # If user leaves pseudopotential empty, assign a default based on element.
        if not pseudo:
            pseudo = f"{elem}.pbesol.UPF"
            color_print(f"No pseudopotential provided for {elem}. Using default: {pseudo}", Fore.YELLOW)
        atomic_species.append({"element": elem, "mass": mass, "pseudopotential": pseudo})
    # Build unique atomic species dictionary (to avoid duplicates)
    unique_species = {}
    for sp in atomic_species:
        elem = sp["element"]
        # If already present, you might average masses or keep first entry; here we keep the first.
        if elem not in unique_species:
            unique_species[elem] = sp
    params["atomic_species"] = list(unique_species.values())
    color_print("Unique atomic species: " + ", ".join(unique_species.keys()), Fore.WHITE)

    try:
        num_atoms = int(input("\nEnter the total number of atoms in the unit cell: ").strip())
    except ValueError:
        color_print("Invalid number for total atoms. Exiting.", Fore.RED)
        sys.exit(1)

    atomic_positions = []
    for i in range(num_atoms):
        color_print(f"\nEnter fractional coordinates for atom {i+1}:", Fore.WHITE)
        elem = input("  Element symbol: ").strip().upper()
        try:
            x = float(input("  Fractional coordinate x: ").strip())
            y = float(input("  Fractional coordinate y: ").strip())
            z = float(input("  Fractional coordinate z: ").strip())
        except ValueError:
            color_print("Invalid coordinate input. Exiting.", Fore.RED)
            sys.exit(1)
        atomic_positions.append({"element": elem, "x": x, "y": y, "z": z})
    params["atomic_positions"] = atomic_positions

    band_gap_input = input("\nEnter the band gap (in eV) [Press Enter for default 1.0 eV]: ").strip()
    if band_gap_input == "":
        params["band_gap"] = 1.0
    else:
        try:
            params["band_gap"] = float(band_gap_input)
        except ValueError:
            color_print("Invalid band gap input. Defaulting to 1.0 eV.", Fore.YELLOW)
            params["band_gap"] = 1.0

    recommended_nbnd = compute_recommended_nbnd(params)
    nbnd_input = input(f"\nEnter the number of bands [Press Enter to use recommended {recommended_nbnd}]: ").strip()
    if nbnd_input == "":
        params["nbnd"] = recommended_nbnd
    else:
        try:
            params["nbnd"] = int(nbnd_input)
        except ValueError:
            color_print("Invalid input for nbnd. Using recommended value.", Fore.YELLOW)
            params["nbnd"] = recommended_nbnd

    # Build AI prompt.
    desc_text = f" The material is described as: {params['material_description']}." if params["material_description"] else ""
    props_text = f" Important properties: {params['important_properties']}." if params["important_properties"] else ""
    prompt = (
        f"As your Expert AI Assistant, I'm simulating a crystal material named '{params['material_name']}'."
        + desc_text
        + props_text
        + f" Current parameters: lattice = {params['lattice_parameter']} Å, ecutwfc = {params['ecutwfc']} Ry, "
        + f"SCF k-points = {params['kpoints_scf']}, NSCF k-points = {params['kpoints_nscf']}, nbnd = {params['nbnd']}, "
        + f"band gap = {params['band_gap']} eV. The atomic species and positions have been defined."
        " Please evaluate these settings for accuracy and efficiency, and provide references in ACS Nano style with DOIs."
    )

    print_divider()
    print_header("Current AI Expert Evaluation")
    color_print(prompt, Fore.WHITE)
    print_divider()

    while True:
        color_print("Gathering AI recommendations. Please wait...\n", Fore.WHITE)
        ai_response = query_perplexity(prompt, max_tokens=1024)
        print_header("AI Expert Recommendations and Literature Citations")
        print_ai_response_with_references(ai_response)
        params["ai_recommendations"] = ai_response

        satisfied = input("\nAre you satisfied with these recommendations? (y/n): ").strip().lower()
        if satisfied.startswith("y") or satisfied == "":
            break
        else:
            color_print("\nLet's refine the parameters further.\n", Fore.WHITE)
            new_lat = input(f"Enter updated lattice parameter (Å) [current: {params['lattice_parameter']}]: ").strip()
            if new_lat:
                try:
                    params["lattice_parameter"] = float(new_lat)
                except ValueError:
                    color_print("Invalid input. Retaining old lattice parameter.", Fore.YELLOW)
            new_cut = input(f"Enter updated ecutwfc (Ry) [current: {params['ecutwfc']}]: ").strip()
            if new_cut:
                try:
                    params["ecutwfc"] = float(new_cut)
                except ValueError:
                    color_print("Invalid input. Retaining old cutoff.", Fore.YELLOW)
            new_bg = input(f"Enter updated band gap (eV) [current: {params['band_gap']}]: ").strip()
            if new_bg:
                try:
                    params["band_gap"] = float(new_bg)
                except ValueError:
                    color_print("Invalid input. Retaining old band gap.", Fore.YELLOW)
            recommended_nbnd = compute_recommended_nbnd(params)
            new_nbnd = input(f"Enter updated nbnd [Press Enter for recommended {recommended_nbnd}]: ").strip()
            if new_nbnd:
                try:
                    params["nbnd"] = int(new_nbnd)
                except ValueError:
                    color_print("Invalid input. Using recommended value.", Fore.YELLOW)
                    params["nbnd"] = recommended_nbnd
            else:
                params["nbnd"] = recommended_nbnd

            prompt = (
                f"As your Expert AI Assistant, I'm re-simulating '{params['material_name']}'."
                + desc_text
                + props_text
                + f" Updated parameters: lattice = {params['lattice_parameter']} Å, ecutwfc = {params['ecutwfc']} Ry, "
                + f"SCF k-points = {params['kpoints_scf']}, NSCF k-points = {params['kpoints_nscf']}, nbnd = {params['nbnd']}, "
                + f"band gap = {params['band_gap']} eV."
                " Please re-evaluate these settings with references in ACS Nano style with DOIs."
            )
            print_divider()
            print_header("Current AI Expert Re-Evaluation")
            color_print(prompt, Fore.WHITE)
            print_divider()

    params["formula_pretty"] = params["material_name"].replace(" ", "_")
    print_final_parameters(params)

    try:
        with open(params_file, "w", encoding="utf-8") as f:
            json.dump(params, f, indent=4)
        color_print(f"\nOptimized simulation parameters have been saved to {params_file}.", Fore.WHITE)
    except Exception as e:
        color_print(f"ERROR: Could not save parameters to file: {e}", Fore.RED)
        sys.exit(1)

    try:
        with open(params_file, "r", encoding="utf-8") as f:
            reloaded = json.load(f)
        color_print("\nSuccessfully re-loaded the newly created JSON. Check below:", Fore.WHITE)
        for k, v in reloaded.items():
            color_print(f"  {k}: {v}", Fore.WHITE)
    except Exception as e:
        color_print(f"ERROR: Could not re-load the JSON file: {e}", Fore.RED)

    color_print("\nYou may now proceed with the subsequent simulation steps.", Fore.WHITE)

if __name__ == "__main__":
    main()
