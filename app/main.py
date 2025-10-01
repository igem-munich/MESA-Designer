from pathlib import Path
import sys
import streamlit as st
from stmol import *
from streamlit_js_eval import streamlit_js_eval
from streamlit_sortables import sort_items
import py3Dmol
from annotated_text import annotated_text
import zipfile
import io
import json
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction.Restriction_Dictionary import rest_dict
import dnachisel
from streamlit_downloader import downloader
from streamlit_scroll_navigation import scroll_navbar

# Add the parent directory of the current file to sys.path
# This allows for importing modules from the 'util' package.
current_dir = Path(__file__).resolve().parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

# Import custom utility functions and data from the 'util' package
from util.antibody_search import search_antibodies
from util import TMD_DATA, CTEV_DATA, NTEV_DATA, TEVP_DATA, PRS_DATA, AIP_DATA, FRET_ICDs, CHAIN_COLORS, SIGNAL_SEQS, TAG_SEQS
from util.pdb_interaction import extract_chains_from_pdb, get_pdb_from_rcsb

# Set Streamlit page configuration (must be called before any other Streamlit command)
st.set_page_config(page_title="MESA-Designer", layout="wide", page_icon="🧬")

# Initialize Streamlit's session state.
# Session state is used to persist data across reruns of the Streamlit app.
state = st.session_state

# Initialize session state variables if they don't already exist.
# These variables store various selections, search results, and generated sequences.
if "sabdab" not in state: # Stores results from antibody database search (e.g., SAbDab)
    state.sabdab = None
if "skempi" not in state: # Stores additional affinity data from search results (e.g., SKEMPI)
    state.skempi = None
if "pdbs" not in state: # Stores file paths of found PDB structures (from an earlier offline implementation)
    state.pdbs = None
if "tmds" not in state: # Dictionary to store transmembrane domain (TMD) sequences per chain
    state.tmds = {}
if "linkers" not in state: # Dictionary to store TMD-binder linker sequences per chain
    state.linkers = {}
if "pdb_fasta" not in state: # Stores FASTA format sequences extracted from selected PDB chains
    state.pdb_fasta = None
if "current_pdb_chains_data" not in state: # Stores extracted data (sequences, IDs) for each chain of the selected PDB
    state.current_pdb_chains_data = {}
if "highlight_selection" not in state: # Dictionary to store residue selections for highlighting in the 3D viewer (per chain)
    state.highlight_selection = {}
if "prev_pdb_selection" not in state: # Stores the previously selected PDB ID to detect changes
    state.prev_pdb_selection = None
if "binder_fasta" not in state: # FASTA format string to keep track of the assembled binder sequence
    state.binder_fasta = ""
if "protease_chain_association" not in state: # Dictionary to track which MESA chains are associated with protease parts (N-term, C-term, or complete)
    state.protease_chain_association = {"n": set(), "c": set(), "complete": set()}
if "cargo_chain_association" not in state: # Set to track which MESA chains the cargo sequence should be attached to
    state.cargo_chain_association = set()
if "aip_chain_association" not in state: # Set to track which MESA chains auto-inhibitory peptide (AIP) sequences should be attached to
    state.aip_chain_association = set()
if "tag_chain_association" not in state: # Dictionary to track which detection tags are added to which MESA chains
    state.tag_chain_association = {}
if "construct_list_formatted" not in state: # Dictionary storing formatted construct data (for annotated_text display and GenBank file generation)
    state.construct_list_formatted = {}
if "construct_list_unformatted" not in state: # Dictionary storing unformatted construct sequences (simply the sequence)
    state.construct_list_unformatted = {}
if "download_data" not in state: # Stores a BytesIO object containing the final ZIP file for download
    state.download_data = None
if "prev_search" not in state: # Stores the previous search query to detect changes in the search field
    state.prev_search = ""
if "optimization_settings" not in state: # Stores settings for sequence optimization (e.g., restriction enzymes, species)
    state.optimization_settings = {}
if "themes" not in state: # Dictionary to manage light/dark mode themes and their properties
    state.themes = {
        "current_theme": "dark",
        "refreshed": True,
        "dark": {
            "theme.base": "dark",
            "theme.backgroundColor": "#0E1117",
            "theme.primaryColor": "#7D2593",
            "theme.secondaryBackgroundColor": "#262730",
            "theme.textColor": "#FAFAFA",
            "button_face": "🌞",
            "button_icon": ":material/light_mode:"
        },
        "light":  {
            "theme.base": "light",
            "theme.backgroundColor": "#FFFFFF",
            "theme.primaryColor": "#7D2593",
            "theme.secondaryBackgroundColor": "#F0F2F6",
            "theme.textColor": "#31333F",
            "button_face": "🌜",
            "button_icon": ":material/dark_mode:"
        }
    }

    # Initialize Streamlit's theme configuration to dark mode upon first load
    for key, value in state.themes["dark"].items():
        if key.startswith("theme"):
            st._config.set_option(key, value)
if "current_pdb" not in state: # Stores the content of the currently selected PDB file as a string
    state.current_pdb = ""
if "chain_sequences" not in state: # Dictionary to store the final assembled sequences for Chain A and Chain B (binder only, not complete mesa chain)
    state.chain_sequences = {"Chain A": "", "Chain B": ""}
if "validation_warnings" not in state: # Dictionary to store validation warnings by component for summary display
    state.validation_warnings = {}

# Define preset restriction sites, including common BioBrick and iGEM standards.
# These are used for sequence optimization to avoid certain restriction enzyme recognition sites.
REST_SITES: dict[str, list[str]] = {"iGEM BioBrick Full + SacI, NcoI, KpnI, HindIII": ["AgeI", "AvrII", "BamHI", "BglII", "BsaI", "EcoRI", "HindIII", "KpnI", "NcoI", "NgoMIV", "NheI", "NotI", "PstI", "PvuII", "SacI", "SapI", "SpeI", "XbaI", "XhoI"], "iGEM BioBrick Full": ["AgeI","AvrII","BamHI","BglII","BsaI","EcoRI","NgoMIV","NheI","NotI","PstI","PvuII","SapI","SpeI","XbaI","XhoI"], "RFC10": ["EcoRI","NotI","PstI","SpeI","XbaI"], "RFC1000": ["BsaI", "SapI"], "RFC12": ["AvrII","EcoRI","NheI","NotI","PstI","PvuII","SapI","SpeI","XbaI","XhoI"]}

# define preset legal amino acids, if others are contained in chain sequences, the sequence optimization will fail and thus not be available
LEGAL_AMINO_ACIDS: set[str] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*"}


# Input Validation Functions
def validate_search_query(query: str) -> tuple[bool, str]:
    """
    Validate antigen search query.

    :param query: The search query string to validate.
    :return: Tuple of (is_valid, error_message). error_message is empty string if valid.
    """
    query = query.strip()

    if not query:
        return False, "Search query cannot be empty"

    if len(query) > 100:
        return False, "Search query too long (max 100 characters)"

    # Allow alphanumeric, spaces, hyphens, and common scientific notation
    import re
    if not re.match(r'^[a-zA-Z0-9\s\-\_\.]+$', query):
        return False, "Search query contains invalid characters. Use only letters, numbers, spaces, hyphens, underscores, and periods"

    return True, ""


def validate_protein_sequence(seq: str, chain_name: str) -> tuple[str, list[str]]:
    """
    Validate protein sequence input for binder chains.

    :param seq: The protein sequence to validate.
    :param chain_name: Name of the chain for error messages.
    :return: Tuple of (cleaned_sequence, warnings_list). Returns cleaned sequence and list of warning messages.
    """
    # Remove whitespace and convert to uppercase
    cleaned: str = seq.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    warnings: list[str] = []

    if not cleaned:
        return "", []  # Allow empty (user might only use one chain)

    # Check for maximum length
    if len(cleaned) > 5000:
        warnings.append(f"{chain_name}: Sequence exceeds maximum length (5000 amino acids)")
        cleaned = cleaned[:5000]

    # Standard amino acids + stop codon
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")

    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"{chain_name}: Contains non-standard amino acid codes: {', '.join(sorted(invalid_chars))}. This may affect sequence optimization")

    # Check for internal stop codons
    if "*" in cleaned:
        stop_positions: list[int] = [i for i, aa in enumerate(cleaned) if aa == "*"]
        warnings.append(f"{chain_name}: Contains internal stop codon(s) at position(s) {', '.join(map(str, [p+1 for p in stop_positions]))}. This will terminate translation prematurely")

    return cleaned, warnings


def validate_linker_pattern(pattern: str) -> tuple[str, list[str]]:
    """
    Validate linker pattern (should be short amino acid sequence).

    :param pattern: The linker pattern to validate.
    :return: Tuple of (cleaned_pattern, warnings_list).
    """
    cleaned: str = pattern.strip().upper().replace(" ", "")
    warnings: list[str] = []

    if not cleaned:
        warnings.append("Linker pattern cannot be empty")
        return "", warnings

    if len(cleaned) > 100:
        warnings.append("Linker pattern too long (max 100 amino acids). Use 'Repeats' to create longer linkers")
        cleaned = cleaned[:100]

    # Check for non-standard amino acids
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"Contains non-standard characters: {', '.join(sorted(invalid_chars))}")

    # Check for stop codons in linkers
    if "*" in cleaned:
        warnings.append("Linker pattern contains stop codon (*). This may cause translation termination")

    return cleaned, warnings


def validate_linker_sequence(seq: str, chain_name: str) -> tuple[str, list[str]]:
    """
    Validate full linker sequence (generated or manually entered).

    :param seq: The linker sequence to validate.
    :param chain_name: Name of the chain for warning messages.
    :return: Tuple of (cleaned_sequence, warnings_list).
    """
    cleaned: str = seq.strip().upper().replace(" ", "")
    warnings: list[str] = []

    if not cleaned:
        warnings.append(f"{chain_name} linker: Sequence cannot be empty")
        return "", warnings

    if len(cleaned) > 1000:
        warnings.append(f"{chain_name} linker: Exceeds maximum length (1000 amino acids)")
        cleaned = cleaned[:1000]

    # Check for non-standard amino acids
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"{chain_name} linker: Contains non-standard amino acids: {', '.join(sorted(invalid_chars))}")

    # Check for stop codons
    if "*" in cleaned:
        warnings.append(f"{chain_name} linker: Contains stop codon(s). This may terminate translation")

    # Warn about unusual linker characteristics
    if len(cleaned) < 5:
        warnings.append(f"{chain_name} linker: Very short (<5 AA) - may not provide sufficient flexibility")

    if len(cleaned) > 200:
        warnings.append(f"{chain_name} linker: Very long (>200 AA) - verify this is intentional")

    # Check for hydrophobic-heavy linkers (potential aggregation)
    if len(cleaned) > 0:
        hydrophobic: int = sum(1 for aa in cleaned if aa in "AILMFVPW")
        hydrophobic_ratio: float = hydrophobic / len(cleaned)
        if hydrophobic_ratio > 0.5:
            warnings.append(f"{chain_name} linker: >50% hydrophobic - may cause aggregation issues")

    return cleaned, warnings


def validate_tmd_sequence(seq: str, chain_name: str) -> tuple[str, list[str]]:
    """
    Validate transmembrane domain sequence.

    :param seq: The TMD sequence to validate.
    :param chain_name: Name of the chain for warning messages.
    :return: Tuple of (cleaned_sequence, warnings_list).
    """
    cleaned: str = seq.strip().upper().replace(" ", "")
    warnings: list[str] = []

    if not cleaned:
        warnings.append(f"{chain_name} TMD: Sequence cannot be empty when custom TMD is enabled")
        return "", warnings

    # Enforce upper length limit of 200
    if len(cleaned) > 200:
        warnings.append(f"{chain_name} TMD: Exceeds maximum length (200 amino acids)")
        cleaned = cleaned[:200]

    # Warn about length (typical TMDs are 15-50 amino acids)
    if len(cleaned) < 15:
        warnings.append(f"{chain_name} TMD: Short sequence (<15 AA) - typical TMDs are 15-50 amino acids")

    if len(cleaned) > 50:
        warnings.append(f"{chain_name} TMD: Long sequence (>50 AA) - typical TMDs are 15-50 amino acids")

    # Check for non-standard amino acids
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"{chain_name} TMD: Contains non-standard amino acids: {', '.join(sorted(invalid_chars))}")

    # Check for stop codons
    if "*" in cleaned:
        warnings.append(f"{chain_name} TMD: Contains stop codon(s). This may terminate translation")

    # TMD quality checks (warnings)
    if len(cleaned) > 0:
        # TMDs should be predominantly hydrophobic
        hydrophobic: int = sum(1 for aa in cleaned if aa in "AILMFVPW")
        hydrophobic_ratio: float = hydrophobic / len(cleaned)

        if hydrophobic_ratio < 0.5:
            warnings.append(f"{chain_name} TMD: Low hydrophobicity ({hydrophobic_ratio:.1%}) - typical TMDs are >50% hydrophobic")

        # Check for charged residues (can disrupt membrane insertion)
        charged: int = sum(1 for aa in cleaned if aa in "DEKR")
        charged_ratio: float = charged / len(cleaned)
        if charged_ratio > 0.15:  # >15% charged residues
            warnings.append(f"{chain_name} TMD: High charged residue content ({charged_ratio:.1%}) - may affect membrane insertion")

    return cleaned, warnings


def validate_custom_icd_sequence(seq: str) -> tuple[str, list[str]]:
    """
    Validate custom intracellular domain sequence.
    Most permissive validation since ICDs vary widely.

    :param seq: The ICD sequence to validate.
    :return: Tuple of (cleaned_sequence, warnings_list).
    """
    cleaned: str = seq.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    warnings: list[str] = []

    if not cleaned:
        return "", []  # Allow empty (user switches back to TEV design)

    if len(cleaned) > 5000:
        warnings.append("Custom ICD exceeds recommended length (5000 amino acids)")
        cleaned = cleaned[:5000]

    # Check for non-standard amino acids
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"Contains non-standard amino acids: {', '.join(sorted(invalid_chars))}. This may affect sequence optimization")

    # Check for internal stop codons
    if "*" in cleaned:
        internal_stops: int = cleaned.count("*")
        warnings.append(f"Contains {internal_stops} internal stop codon(s) - this will terminate translation prematurely")

    return cleaned, warnings


def validate_protease_sequence(seq: str, protease_type: str) -> tuple[str, list[str]]:
    """
    Validate protease sequence (N-term, C-term, or complete).

    :param seq: The sequence to validate.
    :param protease_type: "N-terminal", "C-terminal", or "Complete" for error messages.
    :return: Tuple of (cleaned_sequence, warnings_list).
    """
    cleaned: str = seq.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    warnings: list[str] = []

    if not cleaned:
        warnings.append(f"{protease_type} protease: Sequence cannot be empty when custom protease is enabled")
        return "", warnings

    # Upper length limit of 1000
    if len(cleaned) > 1000:
        warnings.append(f"{protease_type} protease: Exceeds maximum length (1000 amino acids)")
        cleaned = cleaned[:1000]

    # Check for non-standard amino acids
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"{protease_type} protease: Contains non-standard amino acids: {', '.join(sorted(invalid_chars))}")

    # Warn about stop codons
    if "*" in cleaned:
        warnings.append(f"{protease_type} protease: Contains stop codon(s). This may terminate translation prematurely")

    return cleaned, warnings


def validate_prs_sequence(seq: str) -> tuple[str, list[str]]:
    """
    Validate protease recognition sequence.
    PRS are typically 6-10 amino acids.

    :param seq: The PRS sequence to validate.
    :return: Tuple of (cleaned_sequence, warnings_list).
    """
    cleaned: str = seq.strip().upper().replace(" ", "")
    warnings: list[str] = []

    if not cleaned:
        warnings.append("PRS sequence cannot be empty when custom PRS is enabled")
        return "", warnings

    if len(cleaned) > 50:
        warnings.append("PRS exceeds maximum length (50 amino acids)")
        cleaned = cleaned[:50]

    # Check for non-standard amino acids
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"Contains non-standard amino acids: {', '.join(sorted(invalid_chars))}")

    # Check for stop codons
    if "*" in cleaned:
        warnings.append("Contains stop codon(s). This may affect protease recognition")

    # Warn about unusually long PRS
    if len(cleaned) > 12:
        warnings.append("Unusually long PRS - most recognition sequences are 6-10 amino acids")

    return cleaned, warnings


def validate_cargo_sequence(seq: str) -> tuple[str, list[str]]:
    """
    Validate cargo sequence (most flexible - can be transcription factors, enzymes, etc.).

    :param seq: The cargo sequence to validate.
    :return: Tuple of (cleaned_sequence, warnings_list).
    """
    cleaned: str = seq.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    warnings: list[str] = []

    if not cleaned:
        return "", []  # Allow empty (cargo is optional)

    if len(cleaned) > 10000:
        warnings.append("Cargo sequence exceeds maximum length (10,000 amino acids)")
        cleaned = cleaned[:10000]

    # Check for non-standard amino acids
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"Contains non-standard amino acids: {', '.join(sorted(invalid_chars))}. This may affect sequence optimization")

    # Warn about internal stop codons
    if "*" in cleaned:
        internal_stops: int = cleaned.count("*")
        warnings.append(f"Contains {internal_stops} internal stop codon(s) - this may terminate translation prematurely")

    # Warn about large cargo
    if len(cleaned) > 2000:
        warnings.append(f"Large cargo ({len(cleaned)} AA) - verify expression and folding will be successful")

    return cleaned, warnings


def validate_aip_sequence(seq: str) -> tuple[str, list[str]]:
    """
    Validate auto-inhibitory peptide sequence.
    AIPs are typically 10-30 amino acids.

    :param seq: The AIP sequence to validate.
    :return: Tuple of (cleaned_sequence, warnings_list).
    """
    cleaned: str = seq.strip().upper().replace(" ", "")
    warnings: list[str] = []

    if not cleaned:
        warnings.append("AIP sequence cannot be empty when custom AIP is enabled")
        return "", warnings

    if len(cleaned) > 50:
        warnings.append("AIP exceeds recommended length (50 amino acids). Typical AIPs are 10-30 AA")
        cleaned = cleaned[:50]

    # Check for non-standard amino acids (warning, not error)
    VALID_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY*")
    invalid_chars: set[str] = set(cleaned) - VALID_AA

    if invalid_chars:
        warnings.append(f"Contains non-standard amino acids: {', '.join(sorted(invalid_chars))}")

    # Warn about stop codons
    if "*" in cleaned:
        warnings.append("Contains stop codon(s). This may terminate translation prematurely")

    return cleaned, warnings


def display_warnings(warnings: list[str], component_key: str) -> None:
    """
    Display validation warnings and store them in session state for summary.

    :param warnings: List of warning messages to display.
    :param component_key: Unique key to identify the component (e.g., "chain_a_binder").
    :return: None
    """
    if warnings:
        # Store warnings in session state for summary
        state.validation_warnings[component_key] = warnings

        # Display warnings immediately
        for warning in warnings:
            st.warning(warning, icon="⚠️")
    else:
        # Clear warnings if sequence is now valid
        if component_key in state.validation_warnings:
            del state.validation_warnings[component_key]


def update_chain_highlight_selection(chain_id_to_toggle: str, current_pdb_selection: str) -> None:
    """
    Updates the session state's `highlight_selection` when a PDB chain's checkbox is toggled.
    If a chain is selected, its entire residue range is added to `highlight_selection`.
    If deselected, the chain is removed from `highlight_selection`.

    :param chain_id_to_toggle: The ID of the chain being toggled (e.g., 'A', 'B').
    :param current_pdb_selection: The ID of the currently selected PDB.
    :return: None
    """
    # Construct the unique key for the checkbox based on PDB and chain ID.
    checkbox_key: str = f"{current_pdb_selection}_checkbox_chain_{chain_id_to_toggle}"

    if state[checkbox_key] and chain_id_to_toggle not in state.highlight_selection: #If the checkbox is now checked and if the chain is not already in the selection, add it.
            # Select the entire range of residues for the chain (1-based indexing).
            state.highlight_selection[chain_id_to_toggle] = list(range(state.current_pdb_chains_data[chain_id_to_toggle]["start"], state.current_pdb_chains_data[chain_id_to_toggle]["end"] + 1))
    elif chain_id_to_toggle in state.highlight_selection: # If the checkbox is now unchecked and if the chain is in the selection, remove it.
        del state.highlight_selection[chain_id_to_toggle]


def update_chain_highlight_selection_residues(chain_id_to_change: str, current_pdb_selection: str) -> None:
    """
    Updates the session state's `highlight_selection` for a specific chain
    based on the residue range entered by the user in a text input.

    :param chain_id_to_change: The ID of the chain whose residues are being modified.
    :param current_pdb_selection: The ID of the currently selected PDB.
    :return: None
    """
    # Construct the unique key for the residue input field.
    input_key: str = f"{current_pdb_selection}_residue_input_chain_{chain_id_to_change}"

    # Only update if the chain is currently selected (i.e., exists in highlight_selection).
    if chain_id_to_change in state.highlight_selection.keys():
        # Parse the 'start:end' string from the text input.
        select_from_str: str
        select_to_str: str
        select_from_str, select_to_str = state[input_key].split(":")
        # Update the highlight selection with the new range of residues.
        state.highlight_selection[chain_id_to_change] = list(range(int(select_from_str), int(select_to_str) + 1))


def update_split_protease_value() -> None:
    """
    Toggles the `split_protease_toggle_value` in the session state.
    This function is called by the `on_change` event of the `split_protease` toggle
    to update its internal state before the app reruns, allowing the label to switch dynamically.

    :return: None
    """
    state.split_protease_toggle_value = not state.split_protease_toggle_value


def update_linker_text_input(chain_id: str) -> None:
    """
    Updates the linker sequence in the session state based on user input from a text area.
    This allows users to either use the generated linker or enter a custom sequence.

    :param chain_id: The ID of the chain for which the linker sequence is being updated.
    :return: None
    """
    # Update the specific linker entry in the `linkers` dictionary with the capitalized input.
    state.linkers[f"{chain_id}_linker"] = state[f"{chain_id}_linker_sequence"].upper()


def generate_download() -> None:
    """
    Generates a ZIP file containing the designed MESA constructs in GenBank format,
    selected PDB structure, and a summary file. This function is called when the
    download button is clicked.
    All constructs chosen by the user are added to the zip file.
    """
    # Create an in-memory byte buffer to store the ZIP file content.
    zip_buffer: io.BytesIO = io.BytesIO()

    # Open the zip buffer in write mode. Files are stored as text (zipfile.ZIP_DEFLATED for compression).
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
       # Iterate through all created constructs stored in session state.
        for key, construct in state.construct_list_formatted.items():
            # Skip constructs that the user has unchecked for download.
            if not state[f"{key}_checkbox"]:
                continue

            # Initialize variables for GenBank record creation.
            # Extract a clean name for the record from the FASTA-like header.
            record_name: str = "_".join(construct[0].replace(">", "").split(" ")[1:]).strip()
            record_sequence: str = ""
            record_features: list[SeqFeature] = []

            # Process each part of the construct (sequence, annotated part, etc.).
            for part in construct[1:]:
                if part == "": # Skip empty parts.
                    continue

                # Handle annotated parts (tuples: (sequence, name, color)).
                if isinstance(part, tuple):
                    if len(part[0]) > 0: # Ensure the sequence part is not empty.
                        # Define a feature qualifier dictionary with the name and translation (if using sequence optimization).
                        qualifiers: dict[str, str] = {"name": part[1], "translation": part[0]} if state.sequence_optimization_toggle else {"name": part[1]}
                        # Add a sequence feature for annotated parts, calculating start/end in base pairs (x3 for DNA when using sequence optimization).
                        record_features.append(SeqFeature(location=FeatureLocation(start=(len(record_sequence) * (3 if state.sequence_optimization_toggle else 1)),
                                                                                   end=(len(record_sequence) + len(part[0])) * (3 if state.sequence_optimization_toggle else 1)),
                                                          type="CDS", # Coding Sequence type
                                                          qualifiers=qualifiers)) # Store name and amino acid translation.
                        record_sequence += part[0] # Append amino acid sequence.

                else: # Handle non-annotated string parts (e.g., linkers not specifically annotated).
                    if len(part) > 0: # Ensure the string part is not empty.
                        record_sequence += part # Append amino acid sequence.

            # Define a clean file name for the GenBank file.
            file_name: str = f"{key.replace(' ', '_').replace(':', '').replace('-', '_').replace('/', '_')}.gb"

            if state.sequence_optimization_toggle:
                # Perform sequence optimization using dnachisel.
                # Reverse translate the amino acid sequence to DNA.
                sequence = dnachisel.reverse_translate(record_sequence)
                # Define constraints: avoid specified restriction sites and enforce translation.
                cons = [dnachisel.AvoidPattern(pat+"_site") for pat in state.optimization_settings["restriction_enzymes"]]
                cons.append(dnachisel.EnforceTranslation())
                # Create a DNA optimization problem.
                problem = dnachisel.DnaOptimizationProblem(sequence = sequence, constraints = cons, objectives = [dnachisel.CodonOptimize(species=state.optimization_settings["species"], method="match_codon_usage")])
                # Resolve constraints and optimize the DNA sequence.
                problem.resolve_constraints(final_check = True)
                problem.optimize()
                # Get the optimized sequence as a SeqRecord.
                recordseq: SeqRecord = problem.to_record()

                # Create the final SeqRecord for GenBank output.
                record: SeqRecord = SeqRecord(recordseq.seq, # Use the optimized DNA sequence.
                                   id=record_name,
                                   name=record_name,
                                   description=record_name,
                                   features=record_features,
                                   annotations={"molecule_type": "DNA"}) # Specify molecule type as DNA.

            else:
                # Create the final SeqRecord for GenBank output.
                record: SeqRecord = SeqRecord(Seq(record_sequence), # Use the unoptimized amino acid sequence.
                                              id=record_name,
                                              name=record_name,
                                              description=record_name,
                                              features=record_features,
                                              annotations={"molecule_type": "PROTEIN"}) # Specify molecule type as protein.

            # Write the GenBank record to a temporary string buffer.
            tmp_file: io.StringIO = io.StringIO()
            SeqIO.write(record, tmp_file, "genbank")
            # Add the GenBank file to the ZIP archive.
            zf.writestr(file_name, tmp_file.getvalue())

        # Add the selected PDB structure to the ZIP file if chosen by the user.
        if state.download_sel_pdb and state.pdbs and state.pdb_selection:
            # Store the PDB file in a subfolder within the ZIP archive.
            zf.writestr(f"selected_pdb/{state.pdb_selection}.pdb", state.current_pdb)

        # Add additional data summary if selected.
        if state["download_additional"]:
            # Construct a summary string with key design information.
            summary_content: str = f"""
                MESA Design Tool Output Summary
                ------------------------------
                Generated on: {datetime.today().strftime('%Y-%m-%d %H:%M:%S')}
                Selected PDB: {state.prev_pdb_selection if state.prev_pdb_selection else 'N/A'}

                Selected Binder FASTA:
                {state.pdb_fasta if state.pdb_fasta else 'N/A'}

                Linker Information:
                {json.dumps(state.linkers, indent=2)}

                TMD Information:
                {json.dumps(state.tmds, indent=2)}
                """
            # Add the summary to the ZIP archive.
            zf.writestr("mesa_design_summary.txt", summary_content.strip())

            # Write SAbDab search results to a CSV file in the ZIP archive.
            s: io.StringIO = io.StringIO()
            state.sabdab.to_csv(s)
            zf.writestr("sabdab_data.csv", s.getvalue())

    # Reset the buffer's position to the beginning after writing all data.
    zip_buffer.seek(0)

    # Store the complete ZIP file content in the session state for download.
    state.download_data = zip_buffer.getvalue()


# cache for better performance as users typically switch back and forth, thus allowing for previous determinations to be re-used
@st.cache_data(show_spinner=False)
def enable_sequence_optimization(sequences: list[str], legal_amino_acids: set[str]) -> bool:
    """
    This function takes a list of sequences and determines whether sequence optimization is possible.
    :param sequences: A list of sequences.
    :param legal_amino_acids: A set of legal amino acid characters.
    :return: True if sequence optimization is possible, False otherwise.
    """
    # check if any of the sequences contain illegal characters
    for seq in sequences:
        seq_letters: set[str] = set(letter for letter in seq)
        for aa in seq_letters:
            if aa not in legal_amino_acids:
                print(f"Sequence optimization is not possible due to illegal character in sequence: {seq}; aa: {aa}")
                return False

    # return true otherwise
    return True


# cache due to frequent repeats of same sequence
@st.cache_data(show_spinner=False)
def get_unannotated_construct_list(annotated_construct_list: dict[str, list[str | tuple[str, str] | tuple[str, str, str]]]) -> dict[str, str]:
    construct_list_unformatted = {}
    for key, construct in annotated_construct_list.items():
        construct_sequence = ""
        for part in construct[1:]:
            if part == "":  # Skip empty parts.
                continue

            # Handle annotated parts (tuples: (sequence, name, color)).
            if isinstance(part, tuple):
                if len(part[0]) > 0:  # Ensure the sequence part is not empty.
                    construct_sequence += part[0]  # Append amino acid sequence.

            else:  # Handle non-annotated string parts (e.g., linkers not specifically annotated).
                if len(part) > 0:  # Ensure the string part is not empty.
                    construct_sequence += part  # Append amino acid sequence.

        construct_list_unformatted[key] = construct_sequence

    return construct_list_unformatted


# Streamlit app configuration and theme management.
def change_theme() -> None:
    """
    Toggles between light and dark themes for the Streamlit application.
    Updates Streamlit's internal configuration options for theme properties.

    :return: None
    """
    # Determine the previous theme.
    previous_theme: str = state.themes["current_theme"]
    # Select the theme dictionary (light if current is dark, dark if current is light).
    tdict: dict[str, str] = state.themes["light"] if state.themes["current_theme"] == "dark" else state.themes["dark"]
    # Apply the new theme settings to Streamlit's configuration.
    for vkey, vval in tdict.items():
        if vkey.startswith("theme"):
            st._config.set_option(vkey, vval)

    # Mark that the theme has been refreshed and update the current theme in state.
    state.themes["refreshed"] = False
    if previous_theme == "dark":
        state.themes["current_theme"] = "light"
    elif previous_theme == "light":
        state.themes["current_theme"] = "dark"


# cache version of get_pdb_from_rcsb
@st.cache_data(show_spinner="Fetching Structure from RCSB PDB...")
def get_cached_pdb_from_rcsb(pdb_id: str) -> str | None:
    """
    This function is simply a wrapper around the get_pdb_from_rcsb function which provides streamlit caching
    :param pdb_id: the pdb id to search for
    :return: the pdb file's content or None
    """
    return get_pdb_from_rcsb(pdb_id)


# update scroll navigation
@st.cache_data(show_spinner=False)
def update_scroll_navigation(transmembrane_design: bool, split_design: bool, protease_release_design: bool, cargo_release_design: bool, valine_design: bool) -> tuple[dict[str, str], list[str]]:
    """
    Due to infrequent updates this function caches its output, this, however, requires all relevant parameters to be passed in.
    :param transmembrane_design: Intracellular Mesa or transmembrane design.
    :param split_design: Split design or different chain design.
    :param protease_release_design: Should the protease cut itself free?
    :param cargo_release_design: Should the cargo be released or stay attached to the protease?
    :param valine_design: is valine currently selected on the A-Mesa chain strand
    :return: A tuple of page anchors and respective Icons for the sidebar.
    """
    # set anchor ids
    # initialize constant anchors
    anchor_ids = {"binder": "Ligand binding site",
                  "linker": "Outer linker",
                  "tmd": "Transmembrane Domain",
                  "icd": "Intracellular Component",
                  "download": "Downloads"}

    # initialize images
    anchor_icons = [open("resources/imgs/ecd.svg").read(),  # binder image
                    open("resources/imgs/ecd_linker.svg").read(),  # linker image
                    ]

    # conditionally change outer linker to Linker and remove tmd
    if not transmembrane_design:
        anchor_ids["linker"] = "Linker"
        anchor_ids.pop("tmd")

    else: # insert tmd if selected
        anchor_icons.append(open(f"resources/imgs/tmd{'_valine' if valine_design else ''}.svg").read())  # insert correct tmd image between ecd_linker and next component

    # insert correct icd image
    anchor_icons.append(open(f"resources/imgs/{'split_tev' if split_design else 'complete_tev'}_{'no_protease_release' if not protease_release_design else 'protease_release'}_{'cargo_release' if cargo_release_design else 'no_cargo_release'}.svg").read())

    # add download image
    anchor_icons.append(open("resources/imgs/download.svg").read())
    return anchor_ids, anchor_icons

# Anchor IDs and icons for sidebar scroll navigation
anchor_ids, anchor_icons = update_scroll_navigation(transmembrane_design=(("transmembrane_mesa" in state and state.transmembrane_mesa) or "transmembrane_mesa" not in state),
                                                    split_design=(("split_protease_toggle" in state and state.split_protease_toggle) or "split_protease_toggle" not in state),
                                                    protease_release_design=("release_protease_toggle" in state and state.release_protease_toggle),
                                                    cargo_release_design=(("release_cargo_toggle" in state and state.release_cargo_toggle) or "release_cargo_toggle" not in state),
                                                    valine_design=("A_tmd_selection" in state and state.A_tmd_selection == "Valine"))

# TODO: create all possible mesa icd svgs

with st.sidebar:
    scroll_navbar(
        list(anchor_ids.values()),
        anchor_labels=None,  # Use anchor_ids as labels
        anchor_icons=anchor_icons,
        auto_update_anchor=True
    )

# Create two columns for the app title and theme change button.
col1, col2 = st.columns([1, 0.1], vertical_alignment="bottom")
with col1:
    st.title("MESA-Design Tool")

# Determine the icon for the theme change button based on the current theme.
btn_icon: str = state.themes["light"]["button_icon"] if state.themes["current_theme"] == "light" else state.themes["dark"]["button_icon"]
with col2:
    # Display the theme change button, calling `change_theme` on click.
    st.button("", icon=btn_icon, on_click=change_theme)

# Rerun the app if the theme has just been changed to apply the new theme settings.
if state.themes["refreshed"] == False:
    state.themes["refreshed"] = True
    st.rerun()

# Get the current browser window width using JavaScript to adjust PDB display size.
page_width = streamlit_js_eval(js_expressions="window.innerWidth", key="WIDTH", want_output=True)

### Target Search Field ################################################################################################
# Allow the user to enter a custom binder sequence instead of searching.
custom_binder: bool = st.toggle(
    "Custom Binder",
    value=False,
    key="custom_binder_toggle",
    help="Enter custom chain sequences for MESA-Chains. These could result from a BindCraft run."
)

st.header("Custom Binder Sequence" if custom_binder else "Find Binding Candidate", anchor="Ligand binding site")

if not custom_binder:
    # Create columns for the search input field and search button.
    col1, col2 = st.columns([1, 0.1])

    # Create the search input field for antigen search.
    with col1:
        search_field: str = st.text_input(label="Antigen-Search",
                                     key="search_input",
                                     label_visibility="collapsed",
                                     placeholder="Search target antigen"
                                     )

    # Create the search button.
    with col2:
        search_button: bool = st.button("", key="search_button",
                                  icon=":material/search:", width=45)

    # Perform search if the search button is clicked or if the search query has changed.
    if search_button or (search_field and state.prev_search != search_field):
        is_valid, error_msg = validate_search_query(search_field)

        if not is_valid:
            st.error(error_msg)
        else:
            # Display a spinner while searching and call the antibody search function.
            with st.spinner(f"Searching for: **{search_field}**"):
                state.sabdab, state.skempi, state.pdbs, state.search_duration = search_antibodies(search_field)
                state.prev_search = search_field # Update previous search to track changes.

# Display search results if available and custom binder is not active.
if state.sabdab is not None and not state.custom_binder_toggle:
    if len(state.sabdab) > 0:
        st.subheader("Select Binder")
        # Display the number of results and search duration.
        st.text(str(len(state.sabdab)) + " results, search took " + str(round(state.search_duration.total_seconds(), 2)) + " s")
        # Display search results in a Streamlit dataframe, allowing single-row selection.
        selection = st.dataframe(state.sabdab, selection_mode="single-row", on_select="rerun", column_config={
            "date": st.column_config.DateColumn(format="YYYY MMM"),
            "resolution": st.column_config.ProgressColumn(format="%f", max_value=5.0),
            "organism": None, # Hide certain columns for clarity.
            "heavy_species": None,
            "light_species": None,
            "model": None,
            "antigen_chain": None,
            "short_header": None,
        },
                                 hide_index=True,
                                 key="sabdab_dataframe")

#        grid_builder = GridOptionsBuilder.from_dataframe(state.sabdab.drop(columns=["organism", "heavy_species", "light_species", "model", "antigen_chain", "short_header"]))
#        grid_builder.configure_selection(selection_mode="single", use_checkbox=True)
#        grid_builder.configure_auto_height(autoHeight=True)
#        grid_builder.configure_side_bar(filters_panel=True, columns_panel=False)
#
#        highglight_code = f"""
#        function(cellClassParams) {{
#            if (cellClassParams.
#        }}
#        """
#
#        grid_options = grid_builder.build()
#
#        AgGrid(data=state.sabdab.drop(columns=["organism", "heavy_species", "light_species", "model", "antigen_chain", "short_header"]), gridOptions=grid_options, key="sabdab_dataframe_new")
#
#        # hightlight search text in sabdab dataframe
#        components.html(f"""
#        <script>
#            console.log("JS: Testing element selection from main window.");
#
#            // Select the AgGrid container by its class name
#            var agGridContainer = parent.document.querySelector(".st-key-sabdab_dataframe_new");
#
#            if (agGridContainer) {{
#                console.log("JS: Successfully found AgGrid container by class! Element:", agGridContainer);
#                values = agGridContainer.querySelector(".ag-cell-value");
#                console.log(values);
#            }} else {{
#                console.log("JS: AgGrid container ('.st-key-sabdab_dataframe_new') NOT found. This might be a timing issue.");
#            }}
#
#            // Return a value (required by streamlit_js_eval)
#            "Element selection test complete";
#        </script>
#        """)
#        #streamlit_js_eval(js_expressions=js_code, key="aggrid_selection_test")

        try:
            # Get the PDB ID from the selected row in the dataframe.
            state.pdb_selection = state.sabdab.iloc[selection["selection"]["rows"]]["pdb"].to_numpy()[0]
            # Retrieve the PDB file content from RCSB PDB.
            state.current_pdb = get_cached_pdb_from_rcsb(state.pdb_selection)
        except:
            # Handle cases where no row is selected or an error occurs.
            state["pdb_selection"] = 0
            state.pdb_selection = None
            state.current_pdb = ""

    else:
        # Display an informational message if no targets are found.
        st.info("""
        No targets were found
        Please check for synonymous names of your target and check again
        You can try to use [BindCraft](https://colab.research.google.com/github/martinpacesa/BindCraft/blob/main/notebooks/BindCraft.ipynb) to create a de-novo binder and use its sequence as a custom binder
        """)

### Display Found Binder Structures ####################################################################################
# This section displays the 3D structure and allows chain/residue selection if a PDB is selected and not in custom binder mode.
if state.pdbs and state.pdb_selection and not state.custom_binder_toggle:
    st.header("Inspect Structures")
    # Create columns for PDB display: left sidebar (0.2), main viewer (1), right sidebar (0.2).
    pdb_cols = st.columns([0.3, 1, 0.2])

    # Display radio selection and chain/residue selection
    with pdb_cols[0]:
        prev_pdb_selection = None # Temporary variable to store previous PDB selection.
        # Re-extract chains and reset selection if the PDB selection changes.
        if  state.pdb_selection:
            # Extract detailed chain data (sequences, IDs) from the current PDB.
            state.current_pdb_chains_data = extract_chains_from_pdb(file_content=state.current_pdb)

            # Generate colors for each chain for visualization.
            state.chain_colors = {}
            for chain_id in state.current_pdb_chains_data.keys():
                # Assign a color from a predefined list (CHAIN_COLORS) to each chain ID.
                state.chain_colors[chain_id] = list(state.chain_colors.values())
            prev_pdb_selection = state.pdb_selection # Update previous selection after processing.
        elif prev_pdb_selection:
            state.pdb_selection = prev_pdb_selection # Restore previous selection if current is empty.

        # Chain selection UI.
        st.subheader("Select Chains")
        # Create checkboxes for each chain found in the PDB.
        for chain_id in state.current_pdb_chains_data.keys():
            # Create sub-columns for chain checkbox and residue input.
            col_pdb_selection_left, col_pdb_selection_right = st.columns(
                2, vertical_alignment="bottom")
            with col_pdb_selection_left:
                # Display checkbox for chain selection.
                st.checkbox(
                    f"Chain {chain_id}",
                    key=f"{state.pdb_selection}_checkbox_chain_{chain_id}",
                    # Set initial value based on whether the chain is already in highlight_selection.
                    value=(chain_id in state.highlight_selection.keys()),
                    on_change=update_chain_highlight_selection, # Callback function on change.
                    args=(chain_id, state.pdb_selection) # Arguments for the callback.
                )
            with col_pdb_selection_right:
                checkbox_key: str = f"{state.pdb_selection}_checkbox_chain_{chain_id}"
                # Only show residue input if the chain's checkbox is checked.
                if state[checkbox_key]:
                    # Text input for specifying residue range (e.g., "1:100").
                    st.text_input(
                        f"Residues in Chain {chain_id}",
                        key=f"{state.pdb_selection}_residue_input_chain_{chain_id}",
                        # Set initial value to the full range of currently selected residues.
                        value=f"{state.highlight_selection[chain_id][0]}:{state.highlight_selection[chain_id][-1]}",
                        on_change=update_chain_highlight_selection_residues, # Callback function on change.
                        args=(chain_id, state.pdb_selection), # Arguments for the callback.
                        label_visibility="collapsed" # Hide the default label for a cleaner UI.
                    )

    with pdb_cols[2]:
        st.subheader("Display Style")
        # Radio buttons to select the display style for the 3D protein model.
        display_style = st.radio(
            label="Display Style",
            options=["line", "cross", "stick", "sphere", "cartoon"],
            index=4, # 'cartoon' is selected by default.
            label_visibility="collapsed"
        )

        # link to rcsb for more information
        st.info(f"View on [RCSB PDB](https://www.rcsb.org/structure/{state.pdb_selection})")

    with pdb_cols[1]:
        # Display selected residues using py3Dmol visualization library.
        # Create a py3Dmol viewer instance with a dynamic width and fixed height.
        view = py3Dmol.view(width=page_width*.66, height=500)
        # Set the background color of the viewer based on the current theme.
        view.setBackgroundColor(state.themes[state.themes["current_theme"]]["theme.backgroundColor"])

        # Add the protein model to the viewer with the selected display style.
        add_model(view, xyz=state.current_pdb, model_style=display_style)

        # Apply specific styling to the highlighted (selected) residues.
        for chain_id in state.highlight_selection.keys():
            for residue_index in state.highlight_selection[chain_id]:
                # Style selected residues with a specific color and arrows.
                # CHAIN_COLORS provides colors per chain for highlighting.
                view.setStyle({"chain": chain_id, "resi": residue_index}, {display_style: {
                              "color": CHAIN_COLORS[chain_id], "arrows": True}})

        # Add hover functionality to display residue information (chain, residue type, residue number).
        view.setHoverable({}, True,"""
            function(atom,viewer,event,container) {
                if(!atom.label) {
                    atom.label = viewer.addLabel(`Chain ${atom.chain}:${atom.resn}:${atom.resi}`,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                }
            }""","""
            function(atom,viewer) {
                if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                }
            }""")
        view.zoomTo() # Zoom the view to fit the entire model.

        # Render the py3Dmol viewer in a Streamlit container.
        with st.container(border=True, gap="small"):
            showmol(view, height=500, width=page_width)

    # Show current FASTA selection to the user if any chains/residues are highlighted.
    if len(state.highlight_selection) != 0:
        fasta: str = ""
        # Initialize structure for draggable components to link variable chains.
        items: list[dict[str, str | list[str]]] = [
            {"header": "Chain A", "items": []}, # Container for chains assigned to output Chain A
            {"header": "Chain B", "items": []}, # Container for chains assigned to output Chain B
            {"header": "Components", "items": []} # Container for available PDB chains
        ]

        # Iterate through PDB entries and add selected chains/residues to the FASTA string.
        for chain_id in state.highlight_selection.keys():
            selection = state.highlight_selection[chain_id]
            fasta += f">{chain_id}\n"
            # Extract the sequence slice based on selected residues (1-based to 0-based conversion).
            fasta += f"{state.current_pdb_chains_data[chain_id]['sequence'][selection[0] - state.current_pdb_chains_data[chain_id]['start']:selection[-1] - state.current_pdb_chains_data[chain_id]["start"] + 1]}\n"

            items[2]["items"].append(chain_id) # Add the chain to the 'Components' list for sorting.

        # Remove the last newline character from the FASTA string.
        fasta = fasta[:-1]
        state.pdb_fasta = fasta # Store the generated FASTA in session state.

        # Show the FASTA preview in a Streamlit code block.
        st.subheader("FASTA Preview")
        fasta_preview = st.code(
            body=fasta,
            height="content",
            language="text",
            wrap_lines=True
        )

        state.binder_fasta = fasta # Also update the general binder FASTA state.

        # Section for linking binder variable chains using a sortable interface.
        st.subheader("Link Variable Chains", help="Connect variable Chains from a natural antibody in any order you wish. They will be connected via a flexible linker. ")

        # Display sortable containers for arranging chains into Chain A and Chain B.
        sorted_chains: list = sort_items(items, multi_containers=True)

        # Reset chain sequences before re-building them.
        state.chain_sequences = {"Chain A": "", "Chain B": ""}

        # Assemble the sequence for "Chain A" based on user's sorted selection.
        for i, chain_id in enumerate(sorted_chains[0]["items"]):
            selection = state.highlight_selection[chain_id]
            # Append the selected chain's sequence segment to Chain A.
            state.chain_sequences["Chain A"] += state.current_pdb_chains_data[chain_id]["sequence"][selection[0] - state.current_pdb_chains_data[chain_id]["start"]:selection[-1] - state.current_pdb_chains_data[chain_id]["start"] + 1]

            # Add a linker sequence between chains, except after the last one.
            if i < len(sorted_chains[0]["items"]) - 1:
                state.chain_sequences["Chain A"] += ("GGGGS" * 3) # A common flexible linker.

        # Assemble the sequence for "Chain B" based on user's sorted selection.
        for i, chain_id in enumerate(sorted_chains[1]["items"]):
            selection = state.highlight_selection[chain_id]
            # Append the selected chain's sequence segment to Chain B.
            state.chain_sequences["Chain B"] += state.current_pdb_chains_data[chain_id]["sequence"][selection[0] - state.current_pdb_chains_data[chain_id]["start"]:selection[-1] - state.current_pdb_chains_data[chain_id]["start"] + 1]

            # Add a linker sequence between chains, except after the last one.
            if i < len(sorted_chains[1]["items"]) - 1:
                state.chain_sequences["Chain B"] += ("GGGGS" * 5) # A common flexible linker.

# This section is activated if the "Custom Binder" toggle is enabled.
elif state.custom_binder_toggle:
    state.pdb_fasta = "" # Clear PDB-derived FASTA as custom sequences are used.

    # Create two columns for Chain A and Chain B custom sequence input.
    binder_chain_columns = st.columns(2)

    with binder_chain_columns[0]:
        st.subheader("Chain A")
        # Text area for user to input custom Chain A sequence.
        chain_a_text_input: str = st.text_area(
            label="Chain A Sequence",
            value="",
            height="content",
            max_chars=5000
        )

        if chain_a_text_input:
            cleaned_a, warnings_a = validate_protein_sequence(chain_a_text_input, "Chain A")
            display_warnings(warnings_a, "chain_a_binder")
            state.chain_sequences["Chain A"] = cleaned_a
        else:
            state.chain_sequences["Chain A"] = ""

    with binder_chain_columns[1]:
        st.subheader("Chain B")
        # Text area for user to input custom Chain B sequence.
        chain_b_text_input: str = st.text_area(
            label="Chain B Sequence",
            value="",
            height="content",
            max_chars=5000
        )

        if chain_b_text_input:
            cleaned_b, warnings_b = validate_protein_sequence(chain_b_text_input, "Chain B")
            display_warnings(warnings_b, "chain_b_binder")
            state.chain_sequences["Chain B"] = cleaned_b
        else:
            state.chain_sequences["Chain B"] = ""

    # Build a FASTA-like preview string for custom binder sequences.
    if len(chain_a_text_input) > 0:
        state.pdb_fasta += f"> Chain A\n{chain_a_text_input}\n"

    if len(chain_b_text_input) > 0:
        state.pdb_fasta += f"> Chain B\n{chain_b_text_input}\n"

### Build linker between Binder and TMD ################################################################################
# This section is enabled if either Chain A or Chain B has a sequence (either from PDB or custom).
if len(state.chain_sequences["Chain A"]) > 0 or len(state.chain_sequences["Chain B"]) > 0:
    st.divider() # Add a visual divider.
    st.header("Linker Design", anchor=anchor_ids["linker"])

    # Display an error if no chains have been selected/provided.
    if len(state.chain_sequences["Chain A"]) < 1 and len(state.chain_sequences["Chain B"]) < 1:
        st.error("Select chains to create linkers")

    # Create linker design interface for each chain that has a sequence.
    for chain_id in [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]:
        # Initialize linker state if not already present.
        if f"{chain_id}_linker" not in state.linkers:
            state.linkers[f"{chain_id}_linker"] = ""

        st.subheader(f"{chain_id} Linker")

        # Create columns for linker pattern input, repeats, and generate button.
        linker_columns = st.columns(3, vertical_alignment="bottom")
        with linker_columns[0]:
            # Text input for the linker pattern (e.g., "GGGS").
            linker_input: str = st.text_input(
                key=f"{chain_id}_linker_pattern",
                label=f"{chain_id} Linker Pattern",
                value="GGGS",
                max_chars=1000
            )

            if linker_input:
                cleaned_pattern, pattern_warnings = validate_linker_pattern(linker_input)
                display_warnings(pattern_warnings, f"{chain_id}_linker_pattern")

        with linker_columns[1]:
            # Number input for how many times the pattern should be repeated.
            linker_repeats: int = st.number_input(
                key=f"{chain_id}_linker_repeats",
                label="Repeats",
                min_value=1,
                max_value=1000,
                value=10
            )
        with linker_columns[2]:
            # Button to generate the linker sequence from the pattern and repeats.
            linker_pattern_generate: bool = st.button(
                key=f"{chain_id}_linker_generate",
                label="Generate"
            )

            # If generate button is clicked, compose the linker sequence.
            if state[f"{chain_id}_linker_generate"]:
                state.linkers[f"{chain_id}_linker"] = state[f"{chain_id}_linker_pattern"] * state[f"{chain_id}_linker_repeats"]

        # Text input to display the generated linker sequence or allow custom input.
        linker_sequence: str = st.text_input(
            key=f"{chain_id}_linker_sequence",
            label=f"{chain_id} Linker Sequence",
            value=state.linkers[f"{chain_id}_linker"],
            on_change=update_linker_text_input, # Callback to update state if custom text is entered.
            args=(chain_id, ),
            max_chars=1000
        )

        if state.linkers.get(f"{chain_id}_linker"):
            cleaned_linker, linker_warnings = validate_linker_sequence(
                state.linkers[f"{chain_id}_linker"],
                chain_id
            )
            display_warnings(linker_warnings, f"{chain_id}_linker_sequence")

### TMD picker #########################################################################################################
# This section allows selection/design of Transmembrane Domain (TMD) sequences.
if len(state.chain_sequences["Chain A"]) > 0 or len(state.chain_sequences["Chain B"]) > 0:
    st.divider() # Add a visual divider.
    # Toggle to enable/disable transmembrane design.
    st.toggle(
        label="Transmembrane Design",
        value=True,
        key="transmembrane_mesa"
    )

    if state.transmembrane_mesa:
        # Toggle to allow custom TMD sequences.
        custom_tmd: bool = st.toggle(
            label="Custom TMD",
            value=False,
            key="custom_tmd_toggle",
            help="Use a custom TMD sequence"
        )

        st.header("TMD Picker" if not custom_tmd else "Custom TMD", anchor = "Transmembrane Domain")

        # Display informational text if not using a custom TMD.
        if not custom_tmd:
            st.info("You can base your choice of TMDs on this [Paper](https://pubmed.ncbi.nlm.nih.gov/33392392/) and its [supplementary data](https://pmc.ncbi.nlm.nih.gov/articles/PMC7759213/#sup1).")

        st.info("Use [TMDock](https://membranome.org/tmdock) to predict TMD agglomeration and relative positioning.")

        # Create columns for TMD selection/input for each active chain.
        tmd_sequence_cols = st.columns(len([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]))

        for i, chain_id in enumerate([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]):
            with tmd_sequence_cols[i]:
                # Initialize TMD state for the current chain.
                if f"{chain_id}_tmd" not in state.tmds:
                    state.tmds[f"{chain_id}_tmd"] = ""

                st.subheader(f"{chain_id} TMD")

                if not custom_tmd:
                    # Radio buttons for selecting a pre-defined TMD from TMD_DATA.
                    tmd_selection: str = st.radio(
                        label="tmd_selection",
                        options=TMD_DATA.keys(), # Options are keys from the TMD_DATA dictionary.
                        horizontal=True,
                        label_visibility="collapsed",
                        key=f"{chain_id[-1]}_tmd_selection"
                    )
                    # Display the sequence of the selected pre-defined TMD.
                    st.code(
                        TMD_DATA[tmd_selection][1],
                        language="text",
                        height="content",
                        wrap_lines=True
                    )
                else:
                    # Text input for entering a custom TMD sequence.
                    tmd_sequence: str = st.text_input(
                        label=f"Enter {chain_id} TMD",
                        value="",
                        max_chars=1000,
                        key=f"{chain_id}_tmd_sequence"
                    )

                    if custom_tmd and tmd_sequence:
                        cleaned_tmd, tmd_warnings = validate_tmd_sequence(tmd_sequence, chain_id)
                        display_warnings(tmd_warnings, f"{chain_id}_tmd")

            # Store the selected/entered TMD sequence in session state.
            state.tmds[f"{chain_id}_tmd"] = TMD_DATA[state[f"{chain_id[-1]}_tmd_selection"]][1] if not custom_tmd else state[f"{chain_id}_tmd_sequence"]

    else:
        state.tmds.clear() # Clear TMDs if transmembrane design is disabled.

### INTRACELLULAR PART DESIGNER ########################################################################################
# This section handles the design of the intracellular component, including proteases, cargo, and AIPs.
if len(state.chain_sequences["Chain A"]) > 0 or len(state.chain_sequences["Chain B"]) > 0:
    st.divider() # Add a visual divider.
    st.header("Intracellular Component", anchor="Intracellular Component")

    # Toggle to choose between a custom Intracellular Domain (ICD) or TEV protease-based design.
    custom_icd: bool = st.toggle(
        label="Custom ICD",
        value=False,
        key="custom_icd_toggle",
        help="You can design your own custom ICD if you want or continue building with the TEV-Protease. Note: If you would like to use the general MESA ICD framework, but use a different protein / recognition sequence you can still continue to use this tool."
    )

    # Custom ICD design input.
    if custom_icd:
        st.subheader("Enter ICD Sequence")
        custom_icd_sequence: str = st.text_area(
            label="Custom ICD",
            value="",
            height="content",
            label_visibility="collapsed",
            key="custom_icd_sequence"
        )

        if custom_icd_sequence:
            cleaned_icd, icd_warnings = validate_custom_icd_sequence(custom_icd_sequence)
            display_warnings(icd_warnings, "custom_icd")

    # TEV-Protease ICD Design section.
    else:
        st.subheader("Protease Design")
        protease_design_container = st.container(border=True) # Container for protease design options.

        with protease_design_container:
            # Initialize `split_protease_toggle_value` in state if it doesn't exist.
            # This controls the label of the split protease toggle.
            if "split_protease_toggle_value" not in state:
                state.split_protease_toggle_value = True

            # Toggle to choose between "Split Protease" or "Separate Chains" design.
            split_protease: bool = st.toggle(
                label=("Split Protease" if state.split_protease_toggle_value else "Separate Chains"),
                value=state.split_protease_toggle_value,
                key="split_protease_toggle",
                help="Split protease: The TEV-Protease or custom protein is split into two chains. Separate Chains: The TEV-Protease is fully assembled on one chain and the protein recognition sequence and cargo are on the other chain",
                on_change=update_split_protease_value # Callback to update the toggle's value and label.
            )

            if split_protease:
                # Toggle to decide if the protease should be released upon receptor dimerization.
                release_protease: bool = st.toggle(
                    label="Release Protease",
                    help="Move PRS before Protease. This cuts PRS chain upon Receptor Dimerization and releases Protease intracellularly",
                    value=False,
                    key="release_protease_toggle"
                )

                # Toggle to decide if cargo should be separated from protease upon receptor dimerization.
                st.toggle(
                    label="Release Cargo",
                    help="Keep protease and cargo connected or separate them upon Receptor Dimerizaiton",
                    value=True,
                    key="release_cargo_toggle"
                )

            # Toggle to allow custom protease sequences.
            custom_protease: bool = st.toggle(
                label="Custom Protease",
                value=False,
                key="custom_protease_toggle",
                help="Enter your own protease or use TEV-Protease"
            )

            # Logic for split protease design.
            if state.split_protease_toggle_value:
                if not custom_protease:
                    # Columns for N-terminus and C-terminus split protease selection.
                    split_n_col, split_c_col = st.columns(2)

                    with split_n_col:
                        st.markdown("#### N-Terminus")

                        # Radio buttons to pick N-terminal TEV-Protease variants.
                        n_protease_selection: str = st.radio(
                            "Pick N-Terminal TEV-Protease",
                            options=NTEV_DATA.keys(),
                            label_visibility="collapsed",
                            key="n_protease_selection"
                        )

                        # Display the sequence of the selected N-terminal protease.
                        n_protease_sequence = st.code(
                            NTEV_DATA[n_protease_selection][1],
                            language="text",
                            wrap_lines=True,
                            height="content"
                        )

                        st.markdown("##### Append to Chain")
                        # Columns for associating N-terminal protease with specific MESA chains.
                        protease_association_cols_n = st.columns(len([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]))
                        for col, chain_id in zip(protease_association_cols_n, [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]):
                            with col:
                                # Checkbox for each chain to associate with N-term protease.
                                checkbox: bool = st.checkbox(
                                    chain_id,
                                    key=f"{chain_id}_protease_association_n"
                                )

                                # Update session state based on checkbox selection.
                                if checkbox:
                                    state.protease_chain_association["n"].add(chain_id)
                                else:
                                    if chain_id in state.protease_chain_association["n"]:
                                        state.protease_chain_association["n"].remove(chain_id)

                    with split_c_col:
                        st.markdown("#### C-Terminus")

                        # Radio buttons to pick C-terminal TEV-Protease variants.
                        c_protease_selection: str = st.radio(
                            "Pick C-Terminal TEV-Protease",
                            options=CTEV_DATA.keys(),
                            label_visibility="collapsed",
                            key="c_protease_selection"
                        )

                        # Display the sequence of the selected C-terminal protease.
                        c_protease_sequence = st.code(
                            CTEV_DATA[c_protease_selection][1],
                            language="text",
                            wrap_lines=True,
                            height="content"
                        )

                        st.markdown("##### Append to Chain")
                        # Columns for associating C-terminal protease with specific MESA chains.
                        protease_association_cols_c = st.columns(len([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]))
                        for col, chain_id in zip(protease_association_cols_c, [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]):
                            with col:
                                # Checkbox for each chain to associate with C-term protease.
                                checkbox: bool = st.checkbox(
                                    chain_id,
                                    key=f"{chain_id}_protease_association_c"
                                )

                                # Update session state based on checkbox selection.
                                if checkbox:
                                    state.protease_chain_association["c"].add(chain_id)
                                else:
                                    if chain_id in state.protease_chain_association["c"]:
                                        state.protease_chain_association["c"].remove(chain_id)
                else: # Custom split protease input.
                    st.markdown("#### Custom Protease")
                    st.info("You can use [SPELL](https://dokhlab.med.psu.edu/spell/login.php) to guide your splitting process!")

                    # Columns for custom N-terminus and C-terminus input.
                    split_n_col, split_c_col = st.columns(2)

                    with split_n_col:
                        st.markdown("##### N-Terminus")
                        # Text area for custom N-terminal protease sequence.
                        n_protease_sequence_entry: str = st.text_area(
                            label="N-Terminus",
                            value="",
                            height="content",
                            max_chars=5000,
                            placeholder="Enter n-terminal protease sequence",
                            label_visibility="collapsed",
                            key="n_protease_sequence_entry"
                        )

                        if n_protease_sequence_entry:
                            cleaned_n, warnings_n = validate_protease_sequence(n_protease_sequence_entry, "N-terminal")
                            display_warnings(warnings_n, "n_protease")

                    with split_c_col:
                        st.markdown("##### C-Terminus")
                        # Text area for custom C-terminal protease sequence.
                        c_protease_sequence_entry: str = st.text_area(
                            label="C-Terminus",
                            value="",
                            height="content",
                            max_chars=5000,
                            placeholder="Enter c-terminal protease sequence",
                            label_visibility="collapsed",
                            key="c_protease_sequence_entry"
                        )

                        if c_protease_sequence_entry:
                            cleaned_c, warnings_c = validate_protease_sequence(c_protease_sequence_entry, "C-terminal")
                            display_warnings(warnings_c, "c_protease")

            # Logic for complete protease design (not split).
            else:
                if not custom_protease:
                    st.markdown("#### Protease Selection")

                    # Columns for protease selection and sequence display.
                    selection_col, sequence_col = st.columns([0.5, 2], vertical_alignment="bottom")

                    with selection_col:
                        # Radio buttons to select from existing complete proteases (TEVp).
                        complete_protease_selection: str = st.radio(
                            label="Select protease",
                            options=TEVP_DATA.keys(),
                            label_visibility="collapsed",
                            key="complete_protease_selection"
                        )

                    with sequence_col:
                        # Display the sequence of the selected complete protease.
                        st.code(
                            TEVP_DATA[complete_protease_selection][1],
                            language="text",
                            wrap_lines=True,
                            height="content"
                        )

                    st.markdown("##### Append to Chain")
                    # Columns for associating complete protease with specific MESA chains.
                    protease_association_cols = st.columns(len([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]))
                    for col, chain_id in zip(protease_association_cols, [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]):
                        with col:
                            # Checkbox for each chain to associate with the complete protease.
                            checkbox = st.checkbox(
                                chain_id,
                                key=f"{chain_id}_protease_association"
                            )

                            # Update session state based on checkbox selection.
                            if checkbox:
                                state.protease_chain_association["complete"].add(chain_id)
                            else:
                                if chain_id in state.protease_chain_association["complete"]:
                                    state.protease_chain_association["complete"].remove(chain_id)
                else: # Custom complete protease input.
                    st.markdown("#### Protease Sequence")

                    # Text area for custom complete protease sequence.
                    complete_protease_sequence: str = st.text_area(
                        label="Protease Sequence",
                        value="",
                        height="content",
                        max_chars=5000,
                        label_visibility="collapsed",
                        placeholder="Enter protease sequence",
                        key="custom_protease_sequence"
                    )

                    if complete_protease_sequence:
                        cleaned_complete, warnings_complete = validate_protease_sequence(complete_protease_sequence, "Complete")
                        display_warnings(warnings_complete, "complete_protease")

        # Cargo design section.
        st.subheader("Cargo Design")
        cargo_design_container = st.container(border=True) # Container for cargo design options.

        with cargo_design_container:
            # Toggle to allow custom Protease Recognition Sequence (PRS).
            custom_prs: bool = st.toggle(
                label="Custom PRS",
                value=False,
                key="custom_prs_toggle",
                help="Enter a custom PRS or use an existing version",
            )

            # Pick from default TEVp PRS variants.
            if not custom_prs:
                st.markdown("#### TEVp PRS Variant Selection")

                # Columns for PRS selection and sequence display.
                prs_selection_col, prs_seq_col = st.columns([0.5, 2], vertical_alignment="bottom")
                with prs_selection_col:
                    # Radio buttons to select a PRS variant.
                    prs_selection: str = st.radio(
                        label="PRS Selection",
                        options=PRS_DATA.keys(),
                        label_visibility="collapsed",
                        key="prs_selection"
                    )

                with prs_seq_col:
                    # Display the sequence of the selected PRS.
                    prs_sequence_display = st.code(
                        PRS_DATA[prs_selection][1],
                        language="text",
                        wrap_lines=True,
                        height="content"
                    )

            # Enter custom PRS sequence.
            else:
                st.markdown("#### PRS Sequence")
                # Text area for custom PRS sequence.
                prs_sequence: str = st.text_area(
                    label="Custom PRS Sequence",
                    value="",
                    height="content",
                    max_chars=5000,
                    label_visibility="collapsed",
                    key="custom_prs_sequence",
                    placeholder="Enter PRS sequence"
                )

                if prs_sequence:
                    cleaned_prs, prs_warnings = validate_prs_sequence(prs_sequence)
                    display_warnings(prs_warnings, "custom_prs")

            # Cargo sequence input.
            st.markdown("#### Cargo Sequence")
            cargo_sequence: str = st.text_area(
                label="Cargo Sequence",
                value="",
                height="content",
                label_visibility="collapsed",
                placeholder="Enter cargo sequence",
                max_chars=10000,
                key="cargo_sequence"
            )

            if cargo_sequence:
                cleaned_cargo, cargo_warnings = validate_cargo_sequence(cargo_sequence)
                display_warnings(cargo_warnings, "cargo")

            # Associate cargo to specific chains.
            st.markdown("##### Append to Chain")
            # Columns for associating cargo with specific MESA chains.
            cargo_chain_cols = st.columns(len([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]))
            for col, chain_id in zip(cargo_chain_cols, [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]):
                with col:
                    # Checkbox for each chain to associate with cargo.
                    checkbox = st.checkbox(
                        label=chain_id,
                        key=f"{chain_id}_cargo_association"
                    )

                    # Update session state based on checkbox selection.
                    if checkbox:
                        state.cargo_chain_association.add(chain_id)
                    else:
                        if chain_id in state.cargo_chain_association:
                            state.cargo_chain_association.remove(chain_id)

    # Further options section for intracellular design.
    st.subheader("Further Options")
    further_options_container = st.container(border=True) # Container for additional options.

    with further_options_container:
        if not custom_icd: # These options are only available if not using a custom ICD.
            # Toggle to include an Auto-Inhibitory Peptide (AIP).
            aip_toggle: bool = st.toggle(
                label="Include AIP",
                value=False,
                key="include_aip",
                help="Include Auto-Inhibitory Peptide. Can reduce background noise by inhibiting Protease"
            )

            if aip_toggle:
                # Toggle to allow custom AIP sequences.
                custom_aip: bool = st.toggle(
                    label="Custom AIP",
                    value=False,
                    key="custom_aip_toggle",
                    help="Choose from variants of TEVp AIPs or create custom"
                )

                # Pick TEVp AIP variant.
                if not custom_aip:
                    st.markdown("#### AIP Selection")
                    # Columns for AIP selection and sequence display.
                    aip_selection_col, aip_seq_col = st.columns([0.5, 2], vertical_alignment="bottom")
                    with aip_selection_col:
                        # Radio buttons to select an AIP variant.
                        aip_selection: str = st.radio(
                            label="AIP Selection",
                            options=AIP_DATA.keys(),
                            label_visibility="collapsed",
                            key="aip_selection"
                        )

                    with aip_seq_col:
                        # Display the sequence of the selected AIP.
                        st.code(
                            AIP_DATA[aip_selection][1],
                            language="text",
                            wrap_lines=True,
                            height="content"
                        )

                # Custom AIP input.
                else:
                    st.markdown("#### Custom AIP")

                    # Text area for custom AIP sequence.
                    aip_sequence: str = st.text_area(
                        label="AIP Sequence",
                        value="",
                        height="content",
                        max_chars=1000,
                        label_visibility="collapsed",
                        placeholder="Enter AIP sequence",
                        key="custom_aip_sequence"
                    )

                    if aip_sequence:
                        cleaned_aip, aip_warnings = validate_aip_sequence(aip_sequence)
                        display_warnings(aip_warnings, "custom_aip")

                # Associate AIP with chains.
                st.markdown("##### Append to Chain")
                # Columns for associating AIP with specific MESA chains.
                aip_chain_cols = st.columns(len([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]))
                for col, chain_id in zip(aip_chain_cols, [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]):
                    with col:
                        # Checkbox for each chain to associate with AIP.
                        checkbox = st.checkbox(
                            chain_id,
                            key=f"{chain_id}_aip_association"
                        )

                        # Update session state based on checkbox selection.
                        if checkbox:
                            state.aip_chain_association.add(chain_id)
                        else:
                            if chain_id in state.aip_chain_association:
                                state.aip_chain_association.remove(chain_id)

                st.divider() # Visual divider for AIP section.

            else:
                state.aip_chain_association.clear() # Clear AIP associations if AIP is disabled.

        # Toggle to automatically create FRET (Förster Resonance Energy Transfer) chains.
        fret_chains_toggle: bool = st.toggle(
            label="Create FRET sequences",
            value=False,
            key="fret_chains_toggle",
            help="Create FRET chains to test selected binders and TMDs"
        )

        # Toggle to add detection tags.
        tag_toggle: bool = st.toggle(
            label="Detection Tags",
            value=False,
            key="tag_toggle",
            help="Prepend (3x)FLAG Tag to Sequences"
        )

        if tag_toggle:
            # Columns for tag selection for each active chain.
            tag_cols = st.columns(len([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]))

            for i, chain_id in enumerate([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]):
                with tag_cols[i]:
                    st.write(f"{chain_id} Tags")
                    for tag in TAG_SEQS.keys(): # Iterate through available tag sequences.
                        # Checkbox for each tag.
                        tag_checkbox: bool = st.checkbox(
                            label=tag,
                            value=False,
                            key=f"{chain_id}_{tag}_tag"
                        )

                        # Initialize chain entry in `tag_chain_association` if not present.
                        if chain_id not in state.tag_chain_association.keys():
                            state.tag_chain_association[chain_id] = {tag: False for tag in TAG_SEQS.keys()}

                        # Update the tag's status for the current chain.
                        state.tag_chain_association[chain_id][tag] = tag_checkbox

        else:
            state.tag_chain_association = {} # Clear tag associations if tags are disabled.

### DOWNLOAD AND OVERVIEW ##############################################################################################
# This section provides an overview of the designed components and the download functionality.
if len(state.chain_sequences["Chain A"]) > 0 or len(state.chain_sequences["Chain B"]) > 0:
    st.divider() # Add a visual divider.
    st.header("Component Overview")
    overview_container = st.container(border=True) # Container for the overview.

    with overview_container:
        # Display selected tags.
        if state.tag_toggle:
            st.subheader("Chain Tags")
            for chain_id in state.tag_chain_association:
                st.markdown(f"#### {chain_id} Tags")
                for tag in [tag for tag in state.tag_chain_association[chain_id] if state.tag_chain_association[chain_id][tag]]:
                    # Display the sequence of each selected tag.
                    st.code(
                        TAG_SEQS[tag][1],
                        language="text",
                        wrap_lines=True,
                        height="content"
                    )

        # Display signal sequence if transmembrane design is active.
        if state.transmembrane_mesa:
            st.subheader("Signal Sequence")
            st.code(
                SIGNAL_SEQS["CD4"][1], # Use a predefined CD4 signal sequence.
                language="text",
                wrap_lines=True,
                height="content"
            )

        # Binder FASTA display.
        st.subheader("Binder Overview")
        st.code(
            state.pdb_fasta, # Display the assembled binder FASTA.
            language="text",
            wrap_lines=True,
            height="content"
        )

        # Linker combination overview.
        st.divider()
        st.subheader(
            f"Binder Linker{'s' if len([chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]) > 1 else ''}")

        for chain_id in [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]:
            st.markdown(f"#### {chain_id} Linker")
            st.code(
                state.linkers[f"{chain_id}_linker"], # Display the linker sequence for each chain.
                language="text",
                wrap_lines=True,
                height="content"
            )

        # TMD overview if transmembrane design is active.
        if state.transmembrane_mesa:
            st.divider()
            st.subheader("TMD Overview")

            for chain_id in [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]:
                st.markdown(f"#### {chain_id} TMD")
                st.code(
                    state.tmds[f"{chain_id}_tmd"], # Display the TMD sequence for each chain.
                    language="text",
                    wrap_lines=True,
                    height="content"
                )

        # ICD overview.
        st.divider()
        if state.custom_icd_toggle:
            st.subheader("Custom ICD")
            st.code(
                state.custom_icd_sequence, # Display the custom ICD sequence.
                language="text",
                wrap_lines=True,
                height="content"
            )
        else:
            if state.split_protease_toggle:
                st.subheader("Split Protease Overview")
                split_protease_overview_cols = st.columns(2)

                with split_protease_overview_cols[0]:
                    st.markdown("#### N-Terminus")
                    st.code(
                        # Display N-terminal protease sequence (predefined or custom).
                        NTEV_DATA[state.n_protease_selection][1] if not state.custom_protease_toggle else state.n_protease_sequence_entry,
                        language="text",
                        wrap_lines=True,
                        height="content"
                    )

                with split_protease_overview_cols[1]:
                    st.markdown("#### C-Terminus")
                    st.code(
                        # Display C-terminal protease sequence (predefined or custom).
                        CTEV_DATA[state.c_protease_selection][1] if not state.custom_protease_toggle else state.c_protease_sequence_entry,
                        language="text",
                        wrap_lines=True,
                        height="content"
                    )

            else:
                st.subheader("Protease Overview")
                st.markdown("#### Protease Sequence")
                st.code(
                    # Display complete protease sequence (predefined or custom).
                    TEVP_DATA[state.complete_protease_selection][1] if not state.custom_protease_toggle else state.custom_protease_sequence,
                    language="text",
                    wrap_lines=True,
                    height="content"
                )

            if state.include_aip:
                st.markdown(
                    f"#### {'Custom ' if state.custom_aip_toggle else ' '}AIP Sequence")
                st.code(
                    # Display AIP sequence (predefined or custom).
                    AIP_DATA[state.aip_selection][1] if not state.custom_aip_toggle else state.custom_aip_sequence,
                    language="text",
                    wrap_lines=True,
                    height="content"
                )

            st.divider()
            st.subheader("Cargo Overview")
            st.markdown(
                f"#### {'Custom ' if state.custom_prs_toggle else ' '}PRS Sequence")
            st.code(
                # Display PRS sequence (predefined or custom).
                PRS_DATA[state.prs_selection][1] if not state.custom_prs_toggle else state.custom_prs_sequence,
                wrap_lines=True,
                language="text",
                height="content"
            )

            st.markdown("#### Cargo Sequence")
            st.code(
                state.cargo_sequence, # Display the cargo sequence.
                language="text",
                wrap_lines=True,
                height="content"
            )

    if state.validation_warnings:
        st.divider()
        st.header("⚠️ Validation Warnings Summary")

        with st.expander("View All Validation Warnings", expanded=False):
            warning_container = st.container(border=True)

            with warning_container:
                st.markdown("**The following warnings were detected in your design:**")
                st.markdown("---")

                for component_key, warnings in state.validation_warnings.items():
                    component_name = component_key.replace("_", " ").title()

                    st.markdown(f"**{component_name}:**")
                    for warning in warnings:
                        st.warning(warning, icon="⚠️")
                    st.markdown("---")

                st.info("""
                    **Note:** These are warnings, not errors. Your design can still be downloaded.
                    However, sequences with warnings may not be suitable for sequence optimization or may have issues during expression.
                """)

    st.header("Downloads", anchor="Downloads")
    download_container = st.container(border=True) # Container for download options.
    with download_container:
        # Assemble annotated constructs based on all selections.
        # This dictionary will hold the final constructs ready for display and GenBank generation.
        construct_list: dict[str, list[str | tuple[str, str] | tuple[str, str, str]]] = {}

        for chain_id in [chain_id for chain_id in state.chain_sequences if len(state.chain_sequences[chain_id]) > 0]:
            # Start building the current chain's construct parts.
            current_chain: list[str | tuple[str, str] | tuple[str, str, str]] = [
                # Initial methionine if not transmembrane (simulating cytoplasmic start).
                "M" if not state.transmembrane_mesa else "",
                # CD4 Signal Sequence if transmembrane design is active.
                (SIGNAL_SEQS["CD4"][1], "CD4 Signal Sequence", "#74C30EFF") if state.transmembrane_mesa else "",
                # Binder sequence.
                (state.chain_sequences[chain_id], "Binder", "#534cb3"),
                # Linker sequence.
                (state.linkers[f"{chain_id}_linker"], "Linker", "#eba814"),
                # TMD sequence if transmembrane design is active.
                (state.tmds[f"{chain_id}_tmd"], "TMD", "#69ad52") if state.transmembrane_mesa else "",
                # Fixed linker segment.
                "GGGSGGGS" if state.transmembrane_mesa else "",
            ]

            # Prepend tags if enabled.
            if state.tag_toggle:
                for tag in [tag for tag in state.tag_chain_association[chain_id] if state.tag_chain_association[chain_id][tag]]:
                    # Insert tags at position 2 (after initial M/Signal Peptide, before Binder).
                    current_chain.insert(2, (TAG_SEQS[tag][1], f"{tag} Tag", "#26B771FF"))

            # Logic for split protease design.
            if state.split_protease_toggle:
                if chain_id in state.protease_chain_association["n"]:
                    # Construct N-term protease chain.
                    construct_list[f"{chain_id}_N-Term Protease"] = ([f"> {chain_id}_N-Term Protease{'_CARGO' if chain_id in state.cargo_chain_association else ''}\n\n"]
                                                                     + current_chain
                                                                     # PRS if protease release is toggled.
                                                                     + [(PRS_DATA[state.prs_selection][1] if not state.custom_prs_toggle else state.custom_prs_sequence, "PRS", "#b4774b") if state.release_protease_toggle else ""]
                                                                     + ["GGGSGGGS" if state.release_protease_toggle else ""]
                                                                     # N-term protease sequence.
                                                                     + [(NTEV_DATA[state.n_protease_selection][1] if not state.custom_protease_toggle else state.n_protease_sequence_entry, "N-Term Protease", "#bfbd40")]
                                                                     # Linker and PRS/Cargo if release cargo is toggled and cargo is associated.
                                                                     + ["GGGSGGGS" if chain_id in state.cargo_chain_association and state.release_cargo_toggle else ""]
                                                                     + [(PRS_DATA[state.prs_selection][1] if not state.custom_prs_toggle else state.custom_prs_sequence, "PRS", "#b4774b") if chain_id in state.cargo_chain_association and state.release_cargo_toggle else ""]
                                                                     + ["GGGSGGGS" if chain_id in state.cargo_chain_association else ""]
                                                                     + [(state.cargo_sequence, "Cargo", "#bd4258") if chain_id in state.cargo_chain_association else ""]
                                                                     + ["*" if not chain_id in state.aip_chain_association else ""] # Stop codon if no AIP
                                                                     )

                    # Append AIP if associated with this chain.
                    if chain_id in state.aip_chain_association:
                        construct_list[f"{chain_id}_N-Term Protease"].append("GGGSGGGS")
                        construct_list[f"{chain_id}_N-Term Protease"].append((AIP_DATA[state.aip_selection][1] if not state.custom_aip_toggle else state.custom_aip_sequence, "AIP", "#5aa56b"))
                        construct_list[f"{chain_id}_N-Term Protease"].append("*")

                elif chain_id in state.protease_chain_association["c"]:
                    # Construct C-term protease chain.
                    construct_list[f"{chain_id}_C-Term Protease"] = ([f"> {chain_id}_C-Term Protease{'_CARGO' if chain_id in state.cargo_chain_association else ''}\n\n"]
                                                                     + current_chain
                                                                     # PRS if protease release is toggled.
                                                                     + [(PRS_DATA[state.prs_selection][1] if not state.custom_prs_toggle else state.custom_prs_sequence, "PRS", "#b4774b") if state.release_protease_toggle else ""]
                                                                     + ["GGGSGGGS" if state.release_protease_toggle else ""]
                                                                     # C-term protease sequence.
                                                                     + [(CTEV_DATA[state.c_protease_selection][1] if not state.custom_protease_toggle else state.c_protease_sequence_entry, "C-Term Protease", "#3948c6")]
                                                                     # Linker and PRS/Cargo if release cargo is toggled and cargo is associated.
                                                                     + ["GGGSGGGS" if chain_id in state.cargo_chain_association and state.release_cargo_toggle else ""]
                                                                     + [(PRS_DATA[state.prs_selection][1] if not state.custom_prs_toggle else state.custom_prs_sequence, "PRS", "#b4774b") if chain_id in state.cargo_chain_association and state.release_cargo_toggle else ""]
                                                                     + ["GGGSGGGS" if chain_id in state.cargo_chain_association else ""]
                                                                     + [(state.cargo_sequence, "Cargo", "#bd4258") if chain_id in state.cargo_chain_association else ""]
                                                                     + ["*" if not chain_id in state.aip_chain_association else ""] # Stop codon if no AIP
                                                                     )

                    # Append AIP if associated with this chain.
                    if chain_id in state.aip_chain_association:
                        construct_list[f"{chain_id}_C-Term Protease"].append("GGGSGGGS")
                        construct_list[f"{chain_id}_C-Term Protease"].append((AIP_DATA[state.aip_selection][1] if not state.custom_aip_toggle else state.custom_aip_sequence,"AIP", "#5aa56b"))
                        construct_list[f"{chain_id}_C-Term Protease"].append("*")

            else: # Logic for complete protease design.
                if chain_id in state.protease_chain_association["complete"]:
                    # Construct complete protease chain.
                    construct_list[f"{chain_id}_Protease"] = ([f"> {chain_id}_Protease\n\n"]
                                                              + current_chain
                                                              # Complete protease sequence.
                                                              + [(TEVP_DATA[state.complete_protease_selection][1] if not state.custom_protease_toggle else state.custom_protease_sequence, "Protease", "#6b46b9")]
                                                              + ["*" if not chain_id in state.aip_chain_association else ""] # Stop codon if no AIP
                                                              )

                    # Append AIP if associated with this chain.
                    if chain_id in state.aip_chain_association:
                        construct_list[f"{chain_id}_Protease"].append("GGGSGGGS")
                        construct_list[f"{chain_id}_Protease"].append((AIP_DATA[state.aip_selection][1] if not state.custom_aip_toggle else state.custom_aip_sequence, "AIP", "#5aa56b"))
                        construct_list[f"{chain_id}_Protease"].append("*")

                elif chain_id in state.cargo_chain_association:
                    # Construct cargo-only chain (when protease is on another chain).
                    construct_list[f"{chain_id}_Cargo"] = ([f"> {chain_id}_Cargo\n\n"]
                                                           + current_chain
                                                           # PRS, linker, and cargo sequences.
                                                           + [(PRS_DATA[state.prs_selection][1] if not state.custom_prs_toggle else state.custom_prs_sequence, "PRS", "#b4774b"), "GGGSGGGS", (state.cargo_sequence, "Cargo", "#bd4258")]
                                                           + ["*" if not chain_id in state.aip_chain_association else ""] # Stop codon if no AIP
                                                           )

            # Create FRET sequences if enabled.
            if state.fret_chains_toggle:
                construct_list[f"{chain_id}_FRET_mVenus"] = ([f"> {chain_id}_FRET_mVenus\n\n"]
                                                             + current_chain
                                                             # mVenus FRET ICD sequence.
                                                             + [(FRET_ICDs["mVenus"][1], "mVenus", "#43b6bc")]
                                                             + ["*"]
                                                             )

                construct_list[f"{chain_id}_FRET_mCerulean"] = ([f"> {chain_id}_FRET_mCerulean\n\n"]
                                                                + current_chain
                                                                # mCerulean FRET ICD sequence.
                                                                + [(FRET_ICDs["mCerulean"][1], "mCerulean", "#c43b81")]
                                                                + ["*"]
                                                                )

        # Store the assembled formatted construct list in session state.
        state.construct_list_formatted = construct_list

        # Update state with new constructs
        state.construct_list_unformatted = get_unannotated_construct_list(construct_list)

        # Display constructs for user review and selection for download.
        st.subheader("Construct Overview")

        for key, construct in construct_list.items():
            with st.container():
                download_cols = st.columns(
                    [4, 0.1], vertical_alignment="center")
                with download_cols[0]:
                    with st.container(border=True):
                        # Use annotated_text to display the construct with colored annotations.
                        annotated_text(construct)

                with download_cols[1]:
                    # Checkbox to select individual constructs for download.
                    checkbox = st.checkbox(
                        label="Select for download",
                        value=True, # Default to selected.
                        key=f"{key}_checkbox",
                        label_visibility="collapsed"
                    )

        # Sequence Optimization (example: compliance with iGEM regulations)
        st.toggle(
            label="Sequence Optimization",
            value=enable_sequence_optimization(list(state.construct_list_unformatted.values()), LEGAL_AMINO_ACIDS),
            key="sequence_optimization_toggle",
            help="Options to create optimized nucleotide sequences from MESA Chains. Example: Target Organism: Human; Restriction Sites to avoid: BioBrick RFC[10] (iGEM)",
            disabled=not enable_sequence_optimization(list(state.construct_list_unformatted.values()), LEGAL_AMINO_ACIDS),
        )

        if state.sequence_optimization_toggle:
            st.subheader("Sequence Optimization")
            with st.container(border=True):
                optimization_cols = st.columns(2)
                with optimization_cols[0]:
                    # Select box for predefined restriction site avoidance options.
                    option = st.selectbox("Avoid Restriction Sites", ("RFC10", "RFC12", "RFC1000", "iGEM BioBrick Full", "iGEM BioBrick Full + SacI, NcoI, KpnI, HindIII", "Custom"))

                    state.optimization_settings["sel"] = option # Store selected option.
                with optimization_cols[1]:
                    # Select box for species to optimize codon usage for.
                    species = st.selectbox("Select Species to optimize", ("H. Sapiens", "M. Musculus", "D. Melanogaster", "G. Gallus")).replace(". ", "_").lower()
                    state.optimization_settings["species"] = species # Store selected species.

                if state.optimization_settings["sel"] == "Custom":
                    # Multiselect for adding custom restriction enzymes if "Custom" option is chosen.
                    opt = st.multiselect("Add Restriction Enzyme", rest_dict.keys(), key="restriction_add")
                    if opt:
                        state.optimization_settings["custom_enzymes"] = opt # Store custom enzymes.
                        cols = st.columns(3)
                        with cols[0]:
                            st.markdown("#### Restriction Enzyme")
                            for enzyme in state.optimization_settings["custom_enzymes"]:
                                st.write(enzyme) # Display selected enzyme names.
                        with cols[1]:
                            st.markdown("#### Restriction Site")
                            for enzyme in state.optimization_settings["custom_enzymes"]:
                                st.write(rest_dict[enzyme]["site"]) # Display their recognition sites.
                        state.optimization_settings["restriction_enzymes"] = state.optimization_settings["custom_enzymes"]
                else:
                    # If a predefined option is chosen, use its corresponding restriction sites.
                    state.optimization_settings["restriction_enzymes"] = REST_SITES[state.optimization_settings["sel"]]

        # Include additional files in the download.
        st.subheader("Additional Files")
        with st.container(border=True):
            additional_cols = st.columns(2)
            with additional_cols[0]:
                # Checkbox to include PDB structure, disabled if custom binder is used.
                st.checkbox(
                    label="Include PDB Structure",
                    value=not state.custom_binder_toggle,
                    key="download_sel_pdb",
                    disabled=state.custom_binder_toggle
                )

            with additional_cols[1]:
                # Checkbox to include additional data summary and SAbDab results.
                st.checkbox(
                    label="Include additional Data",
                    value=not state.custom_binder_toggle,
                    key="download_additional",
                    disabled=state.custom_binder_toggle
                )

        # Generate download data if constructs are available.

        # Display download button if download data has been generated.
        if st.button("Start Download"):
            generate_download()
            download_complete: bool = downloader(
                data=state.download_data,
                filename="mesa-design.zip",
                content_type="application/zip",
            )

        #if state.download_data:
        #    download: bool = st.download_button(
        #        label="Download Selected",
        #        key="download_selected",
        #        icon=":material/download:",
        #        data=state.download_data, # The BytesIO object containing the ZIP.
        #        file_name="mesa-design.zip",
        #        mime="application/zip"
        #    )
