from pathlib import Path
import json

# Defines the base directory as the parent of the current file.
BASE_DIR: Path = Path(__file__).resolve().parent
# Defines the data directory, located one level up from the base directory in a folder named "data".
DATA_DIR: Path = BASE_DIR.parent / "data"
# Defines the files directory, located one level up from the base directory in a folder named "files".
FILES_DIR: Path = BASE_DIR.parent / "files"
# Defines the resources directory, located one level up from the base directory in a folder named "resources".
RESOURCES_DIR: Path = BASE_DIR.parent / "resources"

# Opens and loads TMD (Transmembrane Domain) data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "tmd/tmd_list.json", "r") as f:
    TMD_DATA: dict[str, list[str]] = dict(json.load(f))

# Opens and loads AIP (Auto-Inhibitory Peptide) data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "aip/aip_list.json", "r") as f:
    AIP_DATA: dict[str, list[str]] = dict(json.load(f))

# Opens and loads FRET (Förster Resonance Energy Transfer) ICDs (Intracellular Domains) data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "FRET/ICDs.json", "r") as f:
    FRET_ICDs: dict[str, list[str]] = dict(json.load(f))

# Opens and loads CTEV (C-Terminal End of TEV Protease) data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "intracellular/CTEV_list.json", "r") as f:
    CTEV_DATA: dict[str, list[str]] = dict(json.load(f))

# Opens and loads NTEV (N-Terminal End of TEV Protease) data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "intracellular/NTEV_list.json", "r") as f:
    NTEV_DATA: dict[str, list[str]] = dict(json.load(f))

# Opens and loads TEVp (TEV protease) data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "intracellular/TEVp_list.json", "r") as f:
    TEVP_DATA: dict[str, list[str]] = dict(json.load(f))

# Opens and loads PRS (Protease Recognition Site) data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "prs/prs_list.json", "r") as f:
    PRS_DATA: dict[str, list[str]] = dict(json.load(f))

# Opens and loads chain color data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "colors/chain_colors.json", "r") as f:
    CHAIN_COLORS: dict[str, list[str]] = dict(json.load(f))

# Opens and loads signal sequences data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "signal_seqs/signal_sequences.json", "r") as f:
    SIGNAL_SEQS: dict[str, list[str]] = dict(json.load(f))

# Opens and loads tag sequences data from a JSON file into a dictionary.
with open(RESOURCES_DIR / "tags/tag_sequences.json", "r") as f:
    TAG_SEQS: dict[str, list[str]] = dict(json.load(f))