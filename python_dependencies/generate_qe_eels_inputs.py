#!/usr/bin/env python3
import os
import re
import sys

def extract_parameters(scf_path):
    """
    Extract the 'prefix' and 'outdir' parameters from the scf.in file.
    If not found, default values are used.
    """
    prefix = None
    outdir = None

    try:
        with open(scf_path, 'r') as f:
            contents = f.read()
    except FileNotFoundError:
        sys.exit(f"Error: The file {scf_path} does not exist.")

    # Look for prefix and outdir in the file (allowing for quotes with either ' or ")
    prefix_match = re.search(r"prefix\s*=\s*['\"]([^'\"]+)['\"]", contents)
    outdir_match = re.search(r"outdir\s*=\s*['\"]([^'\"]+)['\"]", contents)

    if prefix_match:
        prefix = prefix_match.group(1)
    else:
        print("Warning: 'prefix' not found in scf.in. Defaulting to 'qe'.")
        prefix = "qe"

    if outdir_match:
        outdir = outdir_match.group(1)
    else:
        print("Warning: 'outdir' not found in scf.in. Defaulting to './tmp'.")
        outdir = "./tmp"

    return prefix, outdir

def write_file(filepath, content):
    """Writes the provided content to the file at filepath."""
    try:
        with open(filepath, 'w') as f:
            f.write(content)
        print(f"File written successfully: {filepath}")
    except Exception as e:
        sys.exit(f"Error writing file {filepath}: {e}")

def check_for_paw(scf_path):
    """
    Check the ATOMIC_SPECIES block of scf.in for pseudopotential filenames that include "paw".
    Returns True if any are detected, otherwise False.
    """
    try:
        with open(scf_path, 'r') as f:
            contents = f.read()
    except FileNotFoundError:
        sys.exit(f"Error: The file {scf_path} does not exist.")

    lines = contents.splitlines()
    in_species = False
    for line in lines:
        if line.strip().upper() == "ATOMIC_SPECIES":
            in_species = True
            continue
        if in_species:
            # Stop reading if we hit an empty line or another block
            if (line.strip() == "" or
                line.strip().upper().startswith("ATOMIC_POSITIONS") or
                line.strip().startswith("&")):
                break
            tokens = line.split()
            if len(tokens) >= 3:
                pseudopot = tokens[2]
                if "paw" in pseudopot.lower():
                    return True
    return False

def matches_element(element, pseudo_filename):
    """
    Returns True if the filename's first token (split by '.', '_' or '-')
    matches the element symbol (case-insensitive). Example:
      - 'O.pbe...' -> first token is 'O'
      - 'Sr_pbe...' -> first token is 'Sr'
      - 'Co.pbe...' -> first token is 'Co'
    """
    base = re.split(r'[._-]', pseudo_filename, maxsplit=1)[0]
    return base.lower() == element.lower()

def list_nonpaw_pseudos(script_dir):
    """
    Return a list of all pseudopotential filenames in ../pseudo
    that do NOT contain 'paw' in the filename (case-insensitive).
    """
    pseudo_dir = os.path.join(script_dir, "..", "pseudo")
    if not os.path.isdir(pseudo_dir):
        print(f"Warning: The pseudo folder '{pseudo_dir}' does not exist. No non-PAW pseudopotentials found.")
        return []
    all_files = os.listdir(pseudo_dir)
    nonpaw = []
    for f in all_files:
        if os.path.isfile(os.path.join(pseudo_dir, f)):
            if "paw" not in f.lower():
                nonpaw.append(f)
    return nonpaw

def modify_scf_for_eels(scf_path, new_outdir):
    """
    Read the original scf.in, automatically list all non-PAW pseudopotentials
    in ../pseudo that match each element (based on the first token),
    replace any PAW lines with a user-chosen or automatically selected non-PAW pseudo,
    remove any Hubbard corrections, then set outdir = new_outdir.

    Write the modified content to _scf_eels.in in the same folder as scf.in.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    nonpaw_list = list_nonpaw_pseudos(script_dir)
    if not nonpaw_list:
        print("No non-PAW pseudopotentials found in ../pseudo. Exiting.")
        sys.exit(1)

    # Read original scf.in
    try:
        with open(scf_path, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        sys.exit(f"Error: The file {scf_path} does not exist.")

    lines = content.splitlines()

    # Locate ATOMIC_SPECIES block
    species_start = None
    for i, line in enumerate(lines):
        if line.strip().upper() == "ATOMIC_SPECIES":
            species_start = i
            break
    if species_start is None:
        sys.exit("Error: ATOMIC_SPECIES block not found in scf.in.")

    # Collect lines belonging to ATOMIC_SPECIES
    species_block = []
    species_end = None
    for i in range(species_start + 1, len(lines)):
        if (lines[i].strip() == "" or
            lines[i].strip().upper().startswith("ATOMIC_POSITIONS") or
            lines[i].strip().startswith("&")):
            species_end = i
            break
        species_block.append(lines[i])
    if species_end is None:
        species_end = len(lines)

    # Replace any PAW pseudopotentials with non-PAW
    modified_species_block = []
    for line in species_block:
        tokens = line.split()
        if len(tokens) < 3:
            modified_species_block.append(line)
            continue
        element, mass, pseudopot = tokens[0], tokens[1], tokens[2]
        if "paw" in pseudopot.lower():
            print(f"\nDetected PAW pseudopotential for element {element}: {pseudopot}")

            element_candidates = [f for f in nonpaw_list if matches_element(element, f)]
            if not element_candidates:
                sys.exit(f"No non-PAW pseudopotential found for element {element} in ../pseudo. Exiting.")

            if len(element_candidates) == 1:
                chosen_pseudo = element_candidates[0]
                print(f"  Found exactly one match for element {element}: {chosen_pseudo}")
                print("  Using this pseudopotential automatically.")
            else:
                print("Listing available non-PAW pseudopotentials for this element in ../pseudo:")
                for idx, candidate in enumerate(element_candidates, start=1):
                    print(f"   [{idx}] {candidate}")
                choice = input(f"Select a pseudopotential for {element} (1-{len(element_candidates)}): ").strip()
                try:
                    idx_choice = int(choice)
                    if idx_choice < 1 or idx_choice > len(element_candidates):
                        raise ValueError
                except ValueError:
                    sys.exit("Invalid choice. Exiting.")
                chosen_pseudo = element_candidates[idx_choice - 1]
                print(f"Using {chosen_pseudo} for element {element}.")

            new_line = f"{element} {mass} {chosen_pseudo}"
            modified_species_block.append(new_line)
        else:
            modified_species_block.append(line)

    # Rebuild content with updated pseudopotential lines
    new_content_lines = (
        lines[:species_start+1] +
        modified_species_block +
        lines[species_end:]
    )

    # Remove Hubbard corrections
    no_hubbard_lines = []
    found_hubbard = False
    for ln in new_content_lines:
        low = ln.lower()
        if ("lda_plus_u" in low or
            "hubbard_u" in low or
            "hubbard_v" in low or
            "hubbard_j" in low or
            "hubbard_alpha" in low or
            "hubbard_beta" in low):
            found_hubbard = True
            # skip this line
        else:
            no_hubbard_lines.append(ln)

    if found_hubbard:
        print("===========================================================")
        print("Detected Hubbard corrections in scf.in. Removing them so that")
        print("TDDFT/EELS can proceed without the Hubbard error.")
        print("===========================================================")

    # Convert back to string
    new_content = "\n".join(no_hubbard_lines)

    # Force outdir = new_outdir in &CONTROL
    new_content = re.sub(
        r"(outdir\s*=\s*['\"])[^'\"]+(['\"])",
        rf"\1{new_outdir}\2",
        new_content
    )

    # Write the new scf file
    inputs_dir = os.path.dirname(scf_path)
    new_scf_path = os.path.join(inputs_dir, "_scf_eels.in")
    try:
        with open(new_scf_path, 'w') as f:
            f.write(new_content)
    except Exception as e:
        sys.exit(f"Error writing file {new_scf_path}: {e}")

    print(f"\nModified SCF file for EELS generated: {new_scf_path}")
    prefix, _ = extract_parameters(new_scf_path)
    return new_scf_path, prefix

def run_scf_recalculation(scf_input_path, outputs_dir):
    """
    Runs the SCF recalculation using the modified SCF input file,
    automatically using QE from ../qe/pw.exe (no prompt for path).
    """
    np_val = input("Enter number of MPI processes (NP) [default: 4]: ").strip()
    if not np_val:
        np_val = "4"

    script_dir = os.path.dirname(os.path.abspath(__file__))
    qe_exec_dir = os.path.join(script_dir, "..", "qe")
    qe_pw_exe   = os.path.join(qe_exec_dir, "pw.exe")

    os.makedirs(outputs_dir, exist_ok=True)
    output_file = os.path.join(outputs_dir, "called_scf_eels.out")

    cmd = f"mpiexec -np {np_val} \"{qe_pw_exe}\" < \"{scf_input_path}\" > \"{output_file}\" 2>&1"
    print("\n======================================================")
    print("Running SCF recalculation for EELS with the command:")
    print(cmd)
    print("======================================================")

    ret = os.system(cmd)
    if ret != 0:
        sys.exit("Error: SCF recalculation for EELS failed. Please check the output.")
    print("SCF recalculation completed successfully.")
    return

def run_mini_scf(scf_input_path, outputs_dir):
    """
    Runs a stand-alone SCF calculation using scf_input_path,
    automatically picking ../qe/pw.exe. Asks only for # of MPI processes.
    """
    np_val = input("Enter number of MPI processes (cores) for SCF [default: 6]: ").strip()
    if not np_val:
        np_val = "6"

    script_dir = os.path.dirname(os.path.abspath(__file__))
    qe_exec_dir = os.path.join(script_dir, "..", "qe")
    qe_pw_exe   = os.path.join(qe_exec_dir, "pw.exe")

    os.makedirs(outputs_dir, exist_ok=True)
    output_file = os.path.join(outputs_dir, "scf.out")

    cmd = f"mpiexec -np {np_val} \"{qe_pw_exe}\" < \"{scf_input_path}\" > \"{output_file}\" 2>&1"
    print("\n======================================================")
    print("Running a stand-alone SCF calculation with the command:")
    print(cmd)
    print("======================================================")

    ret = os.system(cmd)
    if ret != 0:
        sys.exit("Error: SCF calculation failed. Please check the output.")
    print("SCF calculation completed successfully.")
    return

def main():
    """
    Main driver. Generates EELS input files, checks if scf.in uses PAW pseudos,
    removes any Hubbard corrections, and if needed re-runs SCF with non-PAW pseudos.
    Finally, writes qe_eels_low.in, qe_eels_high.in, qe_eels_low_spectrum.in, and
    qe_eels_high_spectrum.in in ../inputs.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    inputs_dir = os.path.join(script_dir, "..", "inputs")
    outputs_dir = os.path.join(script_dir, "..", "outputs")
    scf_path    = os.path.join(inputs_dir, "scf.in")

    # 1) Extract prefix and outdir from scf.in
    prefix, original_outdir = extract_parameters(scf_path)
    print(f"Extracted prefix: {prefix}, outdir: {original_outdir}")

    # 2) Check for PAW. If found, prompt user to pick from relevant non-PAW pseudos
    if check_for_paw(scf_path):
        print("\n======================================================")
        print("PAW pseudopotentials detected in scf.in.")
        print("EELS calculations cannot use PAW pseudopotentials.")
        print("We will attempt to re-run SCF with non-PAW pseudopotentials.")
        print("======================================================\n")

        new_outdir = "./tmp/tmpo_eels"
        print(f"Using outdir = {new_outdir} for the new SCF calculation.")
        new_scf_path, new_prefix = modify_scf_for_eels(scf_path, new_outdir)

        print(f"\nRecalculating SCF using the modified file:\n  {new_scf_path}")
        print(f"SCF output will be saved to: {os.path.join(outputs_dir, 'called_scf_eels.out')}")
        run_scf_recalculation(new_scf_path, outputs_dir)

        final_scf_input = new_scf_path
        prefix = new_prefix
        outdir = new_outdir
    else:
        print("\nNo PAW pseudopotentials detected in scf.in.")
        print("Using the original SCF input and its results. No SCF recalculation needed for EELS.\n")
        final_scf_input = scf_path
        outdir = original_outdir

    # 3) Generate EELS input files
    print("Generating EELS input files in ../inputs ...")

    qe_eels_low_content = f"""&lr_input
  prefix       = '{prefix}'
  outdir       = '{outdir}'
  restart_step = 200
  restart      = .false.
/
&lr_control
  itermax      = 1000
  q1           = 0.1
  q2           = 0.0
  q3           = 0.0
  approximation= 'TDDFT'    ! includes Hartree+XC
  pseudo_hermitian = .true. ! faster Lanczos
/
"""

    qe_eels_high_content = f"""&lr_input
  prefix       = '{prefix}'
  outdir       = '{outdir}'
  restart_step = 200
  restart      = .false.
/
&lr_control
  itermax      = 500
  q1           = 2.5
  q2           = 0.0
  q3           = 0.0
  approximation= 'TDDFT'
  pseudo_hermitian = .true.
/
"""

    qe_eels_low_spectrum_content = f"""&lr_input
  prefix       = '{prefix}'
  outdir       = '{outdir}'
  eels         = .true.
  itermax0     = 500
  itermax      = 500          ! use extrapolation if you want
  extrapolation= 'osc'
  epsil        = 0.05         ! broadening in Ry (~0.68 eV)
  units        = 1            ! 1 => eV
  start        = 0.0
  end          = 50.0
  increment    = 0.1
  verbosity    = 0
/
"""

    qe_eels_high_spectrum_content = f"""&lr_input
  prefix       = '{prefix}'
  outdir       = '{outdir}'
  eels         = .true.
  itermax0     = 500
  itermax      = 500
  extrapolation= 'osc'
  epsil        = 0.1        ! broadening in Ry ~1.36 eV
  units        = 1          ! eV
  start        = 0.0
  end          = 300.0      ! or 500 eV, or 2000 eV, etc.
  increment    = 0.5
  verbosity    = 0
/
"""

    qe_eels_low_path          = os.path.join(inputs_dir, "qe_eels_low.in")
    qe_eels_high_path         = os.path.join(inputs_dir, "qe_eels_high.in")
    qe_eels_low_spectrum_path = os.path.join(inputs_dir, "qe_eels_low_spectrum.in")
    qe_eels_high_spectrum_path= os.path.join(inputs_dir, "qe_eels_high_spectrum.in")

    write_file(qe_eels_low_path, qe_eels_low_content)
    write_file(qe_eels_high_path, qe_eels_high_content)
    write_file(qe_eels_low_spectrum_path, qe_eels_low_spectrum_content)
    write_file(qe_eels_high_spectrum_path, qe_eels_high_spectrum_content)

    print("\nEELS input files have been generated:")
    print(f"  {qe_eels_low_path}")
    print(f"  {qe_eels_high_path}")
    print(f"  {qe_eels_low_spectrum_path}")
    print(f"  {qe_eels_high_spectrum_path}")
    print("You can now proceed with the EELS calculations.")

if __name__ == "__main__":
    main()
