#!/usr/bin/env python3
"""
Neuromorphic Calculator - AI Material Discovery Assistant

Assists in identifying materials for neuromorphic computing by querying the Materials Project API and Perplexity AI.
Handles incomplete stoichiometry, retrieves material properties, and generates Quantum ESPRESSO scripts if needed.
All AI responses include citations with publication names and DOIs.
"""

import json
import os
import re
import csv
import time
from openai import OpenAI
from mp_api.client import MPRester

# --- Optional Color Support ---
try:
    from colorama import init, Fore, Style
    init(autoreset=True)
except ImportError:
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

def center_text(text, width=80):
    try:
        width = os.get_terminal_size().columns
    except Exception:
        pass
    return text.center(width)

def print_header(header):
    line = "=" * 80
    color_print(center_text(line), Fore.RED)
    color_print(center_text(header), Fore.RED)
    color_print(center_text(line), Fore.RED)

def normalize_formula(formula):
    """Normalize chemical formulas by removing LaTeX-style subscripts, nicknames, and other formatting."""
    formula = re.sub(r"\\_\{(\d+\.\d+|\d+)\}", r"\1", formula)
    formula = re.sub(r"\\_", "", formula)
    formula = formula.replace("\\(", "").replace("\\)", "").replace("\\", "")
    subscripts = {
        "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
        "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9"
    }
    normalized = "".join(subscripts.get(ch, ch) for ch in formula)
    normalized = re.sub(r"\s*\([^)]+\)", "", normalized)
    normalized = normalized.replace("{", "").replace("}", "").replace(" ", "")
    normalized = normalized.replace("(x)", "").replace("(1-x)", "")
    return normalized.strip()

# --- Perplexity API Setup ---
def setup_perplexity_client():
    api_key_file = "./api_keys/api_key_perp.txt"
    if not os.path.isfile(api_key_file):
        color_print(f"Error: {api_key_file} not found. Please create the file with your Perplexity API key.", Fore.RED)
        exit(1)
    with open(api_key_file, "r", encoding="utf-8") as f:
        api_key = f.read().strip()
    return OpenAI(api_key=api_key, base_url="https://api.perplexity.ai")

client = setup_perplexity_client()

def query_perplexity(prompt, max_tokens=1024, temperature=0.7):
    messages = [
        {"role": "system", "content": "Be precise and concise. Provide citations with publication names and DOIs."},
        {"role": "user", "content": prompt}
    ]
    try:
        response = client.chat.completions.create(
            model="sonar-pro",
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        err = str(e)
        if "401" in err or "Authorization Required" in err:
            return "Error: 401 Authorization Required. Please check your Perplexity API key."
        return f"API request failed: {e}"

def expand_formula(raw_formula, application):
    """Expand incomplete stoichiometry using AI or user input."""
    normalized = normalize_formula(raw_formula)
    incomplete_patterns = [
        r"-x\b", r"\(x\s*[≈~]\s*(\d+(?:\.\d+)?)\)", r"\((\d+(?:\.\d+)?)\s*<\s*x\s*<\s*(\d+(?:\.\d+)?)\)",
        r"\(_x\)", r"\(_{1-x}\)"
    ]
    is_incomplete = any(re.search(pattern, normalized, flags=re.IGNORECASE) for pattern in incomplete_patterns)

    if is_incomplete:
        prompt = (
            f"The material formula '{normalized}' has incomplete stoichiometry and is intended for {application}. "
            "Provide a list of specific compounds with defined stoichiometries matching the chemical system of '{normalized}'. "
            "For each compound, include the formula, band gap (in eV), and crystal system. Include citations with publication names and DOIs. "
            "Format each entry as: '<number>. <Formula>' followed by a description on the next line starting with a dash. "
            "Do not use subscripts in chemical formulas; write numbers as regular digits."
        )
        ai_response = query_perplexity(prompt)
        color_print(f"\n[AI Response for known examples of {normalized}]:", Fore.CYAN)
        candidates, _ = parse_candidates_formula_celltype(ai_response)
        if candidates:
            return candidates
        color_print(f"[WARNING] No known examples found by AI for {normalized}.", Fore.YELLOW)
        return [{"formula": normalized, "crystal_system": "", "band_gap": None}]

    if re.search(r"-x\b", normalized, flags=re.IGNORECASE):
        default_value = 0
        user_input = input(f"Detected partial stoichiometry in '{normalized}'. Provide a value for 'x' (default: {default_value}): ").strip()
        x_val = float(user_input) if user_input else default_value
        replacement = f"-{x_val}" if x_val != 0 else ""
        new_formula = re.sub(r"-x\b", replacement, normalized, flags=re.IGNORECASE)
        return [{"formula": new_formula, "crystal_system": "", "band_gap": None}]

    if "(" in normalized and ")" in normalized:
        if "(x)" in normalized and "(1-x)" in normalized:
            parts = re.match(r"\((.*?)\)\(x\)\((.*?)\)\(1-x\)", normalized)
            if parts:
                component1, component2 = parts.group(1), parts.group(2)
                prompt = (
                    f"The material formula '{normalized}' is a pseudobinary alloy for {application}. "
                    "Provide a list of specific compounds with defined stoichiometries matching the chemical system of '{component1}-{component2}'. "
                    "For each compound, include the formula, band gap (in eV), and crystal system. Include citations with publication names and DOIs. "
                    "Format each entry as: '<number>. <Formula>' followed by a description on the next line starting with a dash."
                )
                ai_response = query_perplexity(prompt)
                color_print(f"\n[AI Response for known examples of {normalized}]:", Fore.CYAN)
                candidates, _ = parse_candidates_formula_celltype(ai_response)
                if candidates:
                    return candidates
                return [{"formula": normalized, "crystal_system": "", "band_gap": None}]

        approx_match = re.search(r"\(x\s*[≈~]\s*(\d+(?:\.\d+)?)\)", normalized, flags=re.IGNORECASE)
        if approx_match:
            value = approx_match.group(1)
            base = re.sub(r"\s*\(.*\)", "", normalized)
            candidate = re.sub(r"[xX]", value, base)
            return [{"formula": candidate, "crystal_system": "", "band_gap": None}]
        range_match = re.search(r"\((\d+(?:\.\d+)?)\s*<\s*x\s*<\s*(\d+(?:\.\d+)?)\)", normalized, flags=re.IGNORECASE)
        if range_match:
            lower, upper = range_match.group(1), range_match.group(2)
            base = re.sub(r"\s*\(.*\)", "", normalized)
            return [
                {"formula": re.sub(r"[xX]", lower, base), "crystal_system": "", "band_gap": None},
                {"formula": re.sub(r"[xX]", upper, base), "crystal_system": "", "band_gap": None}
            ]
        base = re.sub(r"\s*\(.*\)", "", normalized)
        return [{"formula": base, "crystal_system": "", "band_gap": None}]

    return [{"formula": normalized, "crystal_system": "", "band_gap": None}]

def parse_candidates_formula_celltype(ai_text):
    """Parse AI response to extract formulas, crystal systems, band gaps, and references."""
    candidates = []
    references = []
    valid_crystal_systems = {
        "cubic", "tetragonal", "orthorhombic", "hexagonal", "trigonal",
        "rhombohedral", "monoclinic", "triclinic", "face-centered cubic", "fcc"
    }

    # Normalize the entire AI response to handle subscripts
    ai_text = normalize_formula(ai_text)
    lines = ai_text.splitlines()
    i = 0
    current_formula = None
    description_lines = []
    ref_section = False

    while i < len(lines):
        line = lines[i].strip()

        # Check if we've reached the references section
        if "references" in line.lower() or "[1]" in line or "[2]" in line:
            ref_section = True
            # Extract references
            ref_match = re.match(r"^\[\d+\]\s*(.*)\s*\(DOI:\s*([^\)]+)\)", line)
            if ref_match:
                paper_name = ref_match.group(1).strip()
                doi = ref_match.group(2).strip()
                references.append({"paper": paper_name, "doi": doi})
            i += 1
            continue

        if ref_section:
            i += 1
            continue

        # Match formula lines (e.g., "1. Ge2Sb2Te5" or just "Ge2Sb2Te5")
        formula_match = re.match(r"^(?:\d+\.\s*)?([A-Za-z0-9]+)(?:$|\s)", line, re.IGNORECASE)
        if formula_match and not line.startswith("-"):
            if current_formula and description_lines:
                # Process the previous candidate
                description = " ".join(description_lines).lower()
                crystal_system_match = re.search(
                    r"(?:crystallizes in a|has a|with a|system is|is a)\s*([a-zA-Z-]+(?:\s*(?:cubic|system|structure|phase|alloy))?)(?:[^A-Za-z]|$)",
                    description
                )
                crystal_system = None
                if crystal_system_match:
                    crystal_system = crystal_system_match.group(1).strip()
                    crystal_system = re.sub(r"\s*(crystalline\s*)?(structure|system|phase|alloy)", "", crystal_system).strip()
                    crystal_system = "face-centered cubic" if "face-centered cubic" in crystal_system else crystal_system
                    cs_clean = crystal_system if crystal_system in valid_crystal_systems else ""
                else:
                    cs_clean = ""
                band_gap_match = re.search(
                    r"band gap(?:\s*(?:ranging\s*from|in\s*the\s*range\s*of))?\s*(?:of|is)?\s*(?:approximately|~|about)?\s*(\d+\.\d+|\d+)(?:\s*(?:to|-)\s*(\d+\.\d+|\d+))?\s*ev",
                    description
                )
                band_gap = None
                if band_gap_match:
                    band_gap_start = float(band_gap_match.group(1))
                    band_gap_end = float(band_gap_match.group(2)) if band_gap_match.group(2) else None
                    band_gap = (band_gap_start + band_gap_end) / 2 if band_gap_end else band_gap_start
                candidates.append({
                    "formula": current_formula,
                    "crystal_system": cs_clean,
                    "band_gap": band_gap
                })
                description_lines = []

            current_formula = formula_match.group(1).strip()
            i += 1
        elif line.startswith("-") and current_formula:
            # Collect description lines
            description_lines.append(line[1:].strip())
            i += 1
        else:
            i += 1

    # Process the last candidate if exists
    if current_formula and description_lines:
        description = " ".join(description_lines).lower()
        crystal_system_match = re.search(
            r"(?:crystallizes in a|has a|with a|system is|is a)\s*([a-zA-Z-]+(?:\s*(?:cubic|system|structure|phase|alloy))?)(?:[^A-Za-z]|$)",
            description
        )
        crystal_system = None
        if crystal_system_match:
            crystal_system = crystal_system_match.group(1).strip()
            crystal_system = re.sub(r"\s*(crystalline\s*)?(structure|system|phase|alloy)", "", crystal_system).strip()
            crystal_system = "face-centered cubic" if "face-centered cubic" in crystal_system else crystal_system
            cs_clean = crystal_system if crystal_system in valid_crystal_systems else ""
        else:
            cs_clean = ""
        band_gap_match = re.search(
            r"band gap(?:\s*(?:ranging\s*from|in\s*the\s*range\s*of))?\s*(?:of|is)?\s*(?:approximately|~|about)?\s*(\d+\.\d+|\d+)(?:\s*(?:to|-)\s*(\d+\.\d+|\d+))?\s*ev",
            description
        )
        band_gap = None
        if band_gap_match:
            band_gap_start = float(band_gap_match.group(1))
            band_gap_end = float(band_gap_match.group(2)) if band_gap_match.group(2) else None
            band_gap = (band_gap_start + band_gap_end) / 2 if band_gap_end else band_gap_start
        candidates.append({
            "formula": current_formula,
            "crystal_system": cs_clean,
            "band_gap": band_gap
        })

    return candidates, references

def save_candidates_txt(candidates, filename="candidates.txt"):
    try:
        with open(filename, "w", encoding="utf-8") as f:
            for cand in candidates:
                f.write(f"{cand['formula']}, {cand['crystal_system']}\n")
        color_print(f"Candidates saved to {filename}.", Fore.GREEN)
    except Exception as e:
        color_print(f"Error saving candidates to {filename}: {e}", Fore.RED)

def convert_doc(doc):
    return doc.dict() if hasattr(doc, "dict") else doc

def save_accumulated_candidates(candidates, filename="results.txt"):
    try:
        with open(filename, "w", encoding="utf-8") as f:
            for cid, cand in candidates.items():
                if isinstance(cand, dict) and cand.get("error"):
                    f.write(f"{cid}: {cand['error']}\n")
                elif isinstance(cand, dict) and cand.get("all_entries"):
                    f.write(f"{cid}: No exact match found. Nearest match:\n")
                    for entry in cand["all_entries"][:1]:
                        entry = convert_doc(entry)
                        f.write(
                            f"  MP-ID: {entry.get('material_id', 'N/A')}, "
                            f"Formula: {entry.get('formula_pretty', 'N/A')}, "
                            f"Band Gap: {entry.get('band_gap', 'N/A')} eV, "
                            f"Crystal System: {entry.get('symmetry', {}).get('crystal_system', 'N/A')}\n"
                        )
                elif isinstance(cand, dict) and cand.get("ai_info"):
                    f.write(f"{cid}: AI info provided:\n  {cand['ai_info']}\n")
                elif isinstance(cand, dict) and cand.get("qe_params"):
                    f.write(f"{cid}: QE parameters generated via AI:\n  {json.dumps(cand['qe_params'], indent=2)}\n")
                else:
                    f.write(
                        f"{cid}: MP-ID {cand.get('material_id', 'N/A')}, "
                        f"Formula: {cand.get('formula_pretty', 'N/A')}, "
                        f"Band Gap: {cand.get('band_gap', 'N/A')} eV, "
                        f"Crystal System: {cand.get('symmetry', {}).get('crystal_system', 'N/A')}\n"
                    )
        color_print(f"Latest search results saved to {filename}.", Fore.GREEN)
    except Exception as e:
        color_print(f"Error saving accumulated candidates to {filename}: {e}", Fore.RED)

class MaterialSearcher:
    BAND_GAP_TOLERANCE = 0.5

    def __init__(self, api_key):
        self.api_key = api_key

    @staticmethod
    def crystal_system_range(system):
        ranges = {
            "triclinic": (1, 2), "monoclinic": (3, 15), "orthorhombic": (16, 74),
            "tetragonal": (75, 142), "hexagonal": (143, 194), "cubic": (195, 230),
            "rhombohedral": (143, 167), "trigonal": (143, 167),
            "face-centered cubic": (195, 230), "fcc": (195, 230)
        }
        return ranges.get(system.lower(), (1, 230))

    def get_all_candidates(self, formula, band_gap=None):
        normalized_formula = normalize_formula(formula)
        fields = ["material_id", "formula_pretty", "band_gap", "symmetry"]
        max_retries = 3
        for attempt in range(max_retries):
            try:
                with MPRester(self.api_key) as mpr:
                    if band_gap is not None:
                        band_gap_range = (max(0, band_gap - self.BAND_GAP_TOLERANCE), band_gap + self.BAND_GAP_TOLERANCE)
                        docs = mpr.materials.summary.search(
                            formula=normalized_formula,
                            band_gap=band_gap_range,
                            fields=fields
                        )
                    else:
                        docs = mpr.materials.summary.search(
                            formula=normalized_formula,
                            fields=fields
                        )
                    if not docs:
                        elements = sorted(set(re.findall(r"[A-Z][a-z]?", normalized_formula)))
                        if elements:
                            chemsys = "-".join(elements)
                            with MPRester(self.api_key) as mpr:
                                if band_gap is not None:
                                    docs = mpr.materials.summary.search(
                                        chemsys=chemsys,
                                        band_gap=band_gap_range,
                                        fields=fields
                                    )
                                else:
                                    docs = mpr.materials.summary.search(
                                        chemsys=chemsys,
                                        fields=fields
                                    )
                    return docs
            except Exception as e:
                if attempt < max_retries - 1:
                    color_print(f"[WARNING] Attempt {attempt + 1} failed for {normalized_formula}: {e}. Retrying...", Fore.YELLOW)
                    time.sleep(2)
                    continue
                color_print(f"[ERROR] Failed to access Materials Project database for {normalized_formula} after {max_retries} attempts: {e}", Fore.RED)
                return None

    def filter_by_properties(self, docs, expected_band_gap, expected_system):
        sg_min, sg_max = self.crystal_system_range(expected_system)
        filtered = []
        for doc in docs:
            doc = convert_doc(doc)
            if doc.get("band_gap") is None:
                continue
            gap_diff = abs(doc.get("band_gap") - expected_band_gap) if expected_band_gap is not None else 0
            if expected_band_gap is not None and gap_diff > self.BAND_GAP_TOLERANCE:
                continue
            symmetry = doc.get("symmetry", {})
            sg_num = symmetry.get("number") if isinstance(symmetry, dict) else None
            if sg_num is not None and sg_min <= sg_num <= sg_max:
                filtered.append((doc, gap_diff))
        return filtered

    def search_candidates(self, candidates, application):
        results = {}
        for candidate in candidates:
            candidate_id = f"{candidate['formula']} ({candidate['crystal_system']})"
            color_print(f"\n[INFO] Searching Materials Project for candidate: {candidate_id}", Fore.WHITE)
            docs = self.get_all_candidates(candidate['formula'], band_gap=candidate.get("band_gap"))
            if docs is None:
                color_print(f"[ERROR] Materials Project search returned None for {candidate_id}.", Fore.RED)
                results[candidate_id] = {"error": "Failed to access Materials Project database."}
                continue
            elif not docs:
                color_print(f"[INFO] No match found for your compound. I will generate a proper file for your candidate based on what's reported in literature.", Fore.YELLOW)
                elements = sorted(set(re.findall(r"[A-Z][a-z]?", candidate['formula'])))
                if elements:
                    chemsys = "-".join(elements)
                    fields = ["material_id", "formula_pretty", "band_gap", "symmetry"]
                    try:
                        with MPRester(self.api_key) as mpr:
                            if candidate.get("band_gap") is not None:
                                band_gap_range = (max(0, candidate['band_gap'] - self.BAND_GAP_TOLERANCE), 
                                                candidate['band_gap'] + self.BAND_GAP_TOLERANCE)
                                substitute_docs = mpr.materials.summary.search(
                                    chemsys=chemsys,
                                    band_gap=band_gap_range,
                                    fields=fields
                                )
                            else:
                                substitute_docs = mpr.materials.summary.search(
                                    chemsys=chemsys,
                                    fields=fields
                                )
                        if substitute_docs:
                            nearest_match = convert_doc(substitute_docs[0])
                            color_print(
                                f"[INFO] Closest substitute found: {nearest_match.get('formula_pretty', 'N/A')} "
                                f"(MP-ID: {nearest_match.get('material_id', 'N/A')})",
                                Fore.YELLOW
                            )
                            choice = input(
                                f"\nMaterial '{candidate_id}' not found. (U)se closest substitute or (G)enerate Quantum ESPRESSO script? [U/G]: "
                            ).strip().lower()
                            if choice.startswith('u'):
                                results[candidate_id] = {"all_entries": [nearest_match]}
                                continue
                    except Exception as e:
                        color_print(f"Error finding substitutes for {candidate_id}: {e}", Fore.RED)
                prompt = (
                    f"No data found in Materials Project for {candidate['formula']} "
                    f"(crystal system: {candidate['crystal_system']}) for {application}. "
                    "Provide a JSON object with: lattice_parameter (Angstrom), ecutwfc (Ry), kpoints_scf, kpoints_nscf, band_gap (eV), "
                    "atomic_species (list of dicts with 'element', 'mass', 'pseudopotential'), "
                    "atomic_positions (list of dicts with 'element', 'x', 'y', 'z'). "
                    "Include references with publication names and DOIs in a 'references' field."
                )
                ai_response = query_perplexity(prompt)
                color_print(f"[AI Response for {candidate_id}]:", Fore.CYAN)
                ai_response_clean = re.sub(r"```.*?```", "", ai_response, flags=re.DOTALL)
                json_start = ai_response_clean.find('{')
                json_end = ai_response_clean.rfind('}')
                if json_start != -1 and json_end != -1:
                    json_str = ai_response_clean[json_start:json_end+1]
                    try:
                        params_obj = json.loads(json_str)
                        results[candidate_id] = {"qe_params": params_obj}
                    except Exception as e:
                        color_print(f"[ERROR] Could not parse JSON: {e}", Fore.RED)
                        results[candidate_id] = {"ai_info": ai_response}
                else:
                    results[candidate_id] = {"ai_info": ai_response}
                continue
            color_print(f"[INFO] Found {len(docs)} entries for {candidate_id}.", Fore.GREEN)
            filtered = self.filter_by_properties(docs, candidate.get("band_gap"), candidate["crystal_system"])
            if filtered:
                best_doc, _ = min(filtered, key=lambda x: x[1])
                results[candidate_id] = best_doc
            else:
                nearest_match = convert_doc(docs[0])
                color_print(
                    f"[INFO] No exact match for {candidate_id}. Closest substitute: {nearest_match.get('formula_pretty', 'N/A')} "
                    f"(MP-ID: {nearest_match.get('material_id', 'N/A')})",
                    Fore.YELLOW
                )
                choice = input(
                    f"\nNo exact match for '{candidate_id}'. (U)se closest substitute or (G)enerate Quantum ESPRESSO script? [U/G]: "
                ).strip().lower()
                if choice.startswith('u'):
                    results[candidate_id] = {"all_entries": [nearest_match]}
                else:
                    prompt = (
                        f"No exact match found in Materials Project for {candidate['formula']} "
                        f"(crystal system: {candidate['crystal_system']}) for {application}. "
                        "Provide a JSON object with: lattice_parameter (Angstrom), ecutwfc (Ry), kpoints_scf, kpoints_nscf, band_gap (eV), "
                        "atomic_species (list of dicts with 'element', 'mass', 'pseudopotential'), "
                        "atomic_positions (list of dicts with 'element', 'x', 'y', 'z'). "
                        "Include references with publication names and DOIs in a 'references' field."
                    )
                    ai_response = query_perplexity(prompt)
                    color_print(f"[AI Response for {candidate_id}]:", Fore.CYAN)
                    ai_response_clean = re.sub(r"```.*?```", "", ai_response, flags=re.DOTALL)
                    json_start = ai_response_clean.find('{')
                    json_end = ai_response_clean.rfind('}')
                    if json_start != -1 and json_end != -1:
                        json_str = ai_response_clean[json_start:json_end+1]
                        try:
                            params_obj = json.loads(json_str)
                            results[candidate_id] = {"qe_params": params_obj}
                        except Exception as e:
                            color_print(f"[ERROR] Could not parse JSON: {e}", Fore.RED)
                            results[candidate_id] = {"ai_info": ai_response}
                    else:
                        results[candidate_id] = {"ai_info": ai_response}
        return results

    @staticmethod
    def load_candidates(filename="candidates.txt"):
        candidates = []
        try:
            with open(filename, "r", newline="", encoding="utf-8") as csvfile:
                reader = csv.reader(csvfile)
                for row in reader:
                    if not row or len(row) < 2:
                        color_print(f"Skipping incomplete row: {row}", Fore.YELLOW)
                        continue
                    formula = row[0].strip()
                    crystal_system = row[1].strip()
                    candidates.append({"formula": formula, "crystal_system": crystal_system, "band_gap": None})
        except Exception as e:
            color_print(f"Error loading candidates from {filename}: {e}", Fore.RED)
        return candidates

def print_candidates(accumulated_candidates, references):
    header = f"{'No.':<4} {'Candidate':<30} {'Material ID':<20} {'Formula':<20} {'Band Gap (eV)':<15} {'Crystal System':<15}"
    color_print("\nPossible Candidate Materials (accumulated):", Fore.WHITE)
    color_print(header, Fore.WHITE)
    color_print("-" * len(header), Fore.WHITE)
    for idx, (cid, cand) in enumerate(accumulated_candidates.items(), start=1):
        if isinstance(cand, dict) and cand.get("error"):
            mat_id = "No MP data"
            formula = "N/A"
            band_gap = "N/A"
            crystal_system = "N/A"
        elif isinstance(cand, dict) and cand.get("all_entries"):
            mat_id = "Nearest match"
            formula = cand["all_entries"][0].get('formula_pretty', cid.split(' (')[0])
            band_gap = str(cand["all_entries"][0].get('band_gap', 'N/A'))
            crystal_system = cand["all_entries"][0].get('symmetry', {}).get('crystal_system', 'N/A')
        elif isinstance(cand, dict) and (cand.get("ai_info") or cand.get("qe_params")):
            mat_id = "AI info"
            formula = cid.split(' (')[0]
            band_gap = "N/A"
            crystal_system = cid.split('(')[1].rstrip(')')
        else:
            mat_id = cand.get('material_id', 'N/A')
            formula = cand.get('formula_pretty', 'N/A')
            band_gap = str(cand.get('band_gap', 'N/A'))
            crystal_system = cand.get('symmetry', {}).get('crystal_system', 'N/A')
        color_print(f"{idx:<4} {cid:<30} {mat_id:<20} {formula:<20} {band_gap:<15} {crystal_system:<15}", Fore.WHITE)
    
    # Print references
    if references:
        color_print("\nReferences:", Fore.CYAN)
        for idx, ref in enumerate(references, start=1):
            color_print(f"{idx}. {ref['paper']}, DOI: {ref['doi']}", Fore.CYAN)

def test_mp_api_key(api_key):
    try:
        with MPRester(api_key) as mpr:
            mpr.materials.summary.search(formula="Si", fields=["material_id"])
        return True
    except Exception as e:
        color_print(f"[ERROR] Materials Project API key validation failed: {e}", Fore.RED)
        return False

def main():
    for fname in ["candidates.txt", "results.txt", "selected_materials.json"]:
        if os.path.exists(fname):
            os.remove(fname)
            color_print(f"[INFO] Cleared previous session file: {fname}", Fore.YELLOW)

    os.system('cls' if os.name == 'nt' else 'clear')
    print_header("Neuromorphic Calculator - AI Material Discovery Assistant")
    color_print("Hi, I'm your AI assistant for neuromorphic material discovery.\n", Fore.WHITE)

    application = input(
        "Describe the neuromorphic application or material keyword:\n"
        "\nExamples:\n"
	"\n"
        " 1. Phase-change material for switching\n"
        " 2. Memristor-based synapse\n"
        " 3. 2D material with STDP\n"
        " 4. Ferroelectric tunnel junction for memory\n"
        " 5. Spintronic device for spike encoding\n"
        " 6. Organic semiconductor for flexible neuromorphic systems\n"
        " 7. Oxide-based resistive switching synapse\n"
        " 8. Superconducting nanowire for low-power computation\n"
        " 9. Magnetic skyrmion for neuron emulation\n"
        "10. Quantum dot array for synaptic plasticity\n"
        "11. Van der Waals heterostructure for neuromorphic interconnects\n"
        "12. Mott insulator for threshold switching\n"
        "13. Carbon nanotube network for neuromorphic computing\n"
        "14. Flexible polymer electrolyte transistor for artificial neurons\n"
        "15. Graphene memristive device for analog memory\n"
	"\n"
        "Your input: "
    ).strip()

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    mp_key_file = os.path.join(project_root, "api_keys", "api_key_mp.txt")
    if not os.path.isfile(mp_key_file):
        color_print(f"Materials Project API key file not found at {mp_key_file}. Exiting.", Fore.RED)
        exit(1)
    with open(mp_key_file, "r", encoding="utf-8") as f:
        mp_api_key = f.read().strip()
    if not test_mp_api_key(mp_api_key):
        exit(1)

    accumulated_candidates = {}
    additional_context = ""
    all_references = []

    while True:
        prompt = f"I need a material for {application}."
        if additional_context:
            prompt += f" The material should also have: {additional_context}."
        prompt += (
            " Provide a numbered list of candidate materials with brief explanations, including the formula, band gap (in eV), and crystal system. "
            "Format each entry as: '<number>. <Formula>' followed by a description on the next line starting with a dash. "
            "Ensure every reference includes the publication name and DOI. Do not use subscripts in chemical formulas."
        )

        color_print("\n[INFO] Searching for candidates...", Fore.GREEN)
        ai_response = query_perplexity(prompt)
        color_print("\n[AI Suggestions]:", Fore.CYAN)
        for line in ai_response.splitlines():
            if any(keyword in line.lower() for keyword in ["doi", "references"]):
                color_print(line, Fore.CYAN)
            else:
                color_print(line, Fore.GREEN)

        simple_candidates = []
        parsed_candidates, references = parse_candidates_formula_celltype(ai_response)
        all_references.extend(references)  # Accumulate references
        for cand in parsed_candidates:
            expanded = expand_formula(cand['formula'], application)
            for exp_cand in expanded:
                exp_cand['crystal_system'] = cand['crystal_system'] if cand['crystal_system'] else exp_cand['crystal_system']
                exp_cand['band_gap'] = cand['band_gap'] if cand['band_gap'] is not None else exp_cand['band_gap']
                simple_candidates.append(exp_cand)

        # Fallback if no candidates are parsed
        if not simple_candidates:
            color_print("[WARNING] No candidate details found in AI response.", Fore.YELLOW)
            manual_input = input("\nWould you like to manually input a candidate formula to search in Materials Project? (Y/N): ").strip().lower()
            if manual_input.startswith('y'):
                formula = input("Enter the chemical formula (e.g., Ge2Sb2Te5): ").strip()
                crystal_system = input("Enter the crystal system (e.g., rhombohedral, cubic): ").strip()
                band_gap_input = input("Enter the band gap in eV (or press Enter to skip): ").strip()
                band_gap = float(band_gap_input) if band_gap_input else None
                simple_candidates.append({
                    "formula": formula,
                    "crystal_system": crystal_system if crystal_system else "",
                    "band_gap": band_gap
                })
            else:
                action = input("\nDo you want to (S)earch for more materials or (E)xit? [S/E]: ").strip().lower()
                if action.startswith('e'):
                    break
                additional_context = input("\nEnter additional characteristics to refine your search: ").strip()
                continue

        save_candidates_txt(simple_candidates)
        candidates_for_search = MaterialSearcher.load_candidates()
        if not candidates_for_search:
            color_print("[WARNING] No valid candidates found in candidates.txt.", Fore.YELLOW)
            action = input("\nDo you want to (S)earch for more materials or (E)xit? [S/E]: ").strip().lower()
            if action.startswith('e'):
                break
            additional_context = input("\nEnter additional characteristics to refine your search: ").strip()
            continue

        searcher = MaterialSearcher(mp_api_key)
        results = searcher.search_candidates(candidates_for_search, application)
        accumulated_candidates.update({k: v for k, v in results.items() if k not in accumulated_candidates})
        save_accumulated_candidates(accumulated_candidates)
        print_candidates(accumulated_candidates, all_references)

        action = input("\nDo you want to (S)earch for more materials, (C)hoose one to simulate, or (E)xit? [S/C/E]: ").strip().lower()
        if action.startswith("c"):
            if not accumulated_candidates:
                color_print("[WARNING] No candidates available. Please search again.", Fore.YELLOW)
                continue
            color_print("\nSelect a candidate by number:", Fore.GREEN)
            candidate_list = list(accumulated_candidates.items())
            for idx, (cid, cand) in enumerate(candidate_list, start=1):
                if isinstance(cand, dict) and cand.get("error"):
                    color_print(f"  {idx}. {cid}: {cand['error']}", Fore.WHITE)
                elif isinstance(cand, dict) and cand.get("all_entries"):
                    color_print(f"  {idx}. {cid}: Nearest match - Formula: {cand['all_entries'][0].get('formula_pretty', 'N/A')}", Fore.WHITE)
                elif isinstance(cand, dict) and (cand.get("ai_info") or cand.get("qe_params")):
                    color_print(f"  {idx}. {cid}: AI provided simulation parameters.", Fore.WHITE)
                else:
                    color_print(f"  {idx}. {cid}: MP-ID {cand.get('material_id', 'N/A')}, "
                                f"Formula: {cand.get('formula_pretty', 'N/A')}", Fore.WHITE)
            try:
                choice = int(input(f"Your choice (1 to {len(candidate_list)}): ").strip())
                if 1 <= choice <= len(candidate_list):
                    selected_candidate_id, selected_candidate = candidate_list[choice - 1]
                    mat_id = selected_candidate.get("material_id", "").strip()
                    if not mat_id or mat_id == "N/A" or "ai_info" in selected_candidate or "qe_params" in selected_candidate:
                        color_print(f"No valid Materials Project ID found for {selected_candidate_id}.", Fore.YELLOW)
                        if "qe_params" in selected_candidate:
                            color_print(f"Using AI-generated simulation parameters for {selected_candidate_id}.", Fore.YELLOW)
                            inputs_folder = os.path.join(project_root, "inputs")
                            os.makedirs(inputs_folder, exist_ok=True)
                            json_file = os.path.join(inputs_folder, "optimized_simulation_parameters.json")
                            with open(json_file, "w", encoding="utf-8") as f:
                                json.dump(selected_candidate["qe_params"], f, indent=2)
                            color_print(f"Simulation parameters saved to {json_file}.", Fore.GREEN)
                            try:
                                import generate_qe_inputs_ad_hoc as qe_ad_hoc
                                qe_ad_hoc.main()
                                color_print(f"Quantum ESPRESSO input files generated successfully for {selected_candidate_id}.", Fore.GREEN)
                            except Exception as e:
                                color_print(f"[ERROR] Failed to generate Quantum ESPRESSO input files: {e}", Fore.RED)
                                continue
                        else:
                            color_print("No valid simulation parameters available. Cannot generate QE inputs.", Fore.RED)
                            continue
                    else:
                        output_file = os.path.join(os.getcwd(), "selected_materials.json")
                        with open(output_file, "w", encoding="utf-8") as f:
                            json.dump({"material_id": mat_id}, f, indent=4)
                        color_print(f"\nSelected candidate saved to {output_file}.", Fore.GREEN)
                        color_print(f"Material ID {mat_id} can be used for further simulation or analysis.", Fore.GREEN)
                    break
                else:
                    color_print(f"Invalid input. Enter a number between 1 and {len(candidate_list)}.", Fore.YELLOW)
            except ValueError:
                color_print("Invalid input. Please enter a valid number.", Fore.YELLOW)
        elif action.startswith("e"):
            color_print("Exiting the Neuromorphic Calculator. Goodbye!", Fore.GREEN)
            break
        else:
            additional_context = input("\nEnter additional characteristics to refine your search: ").strip()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        color_print(f"ERROR: The Neuromorphic Calculator encountered a fatal error: {e}", Fore.RED)
        exit(1)