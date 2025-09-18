import json
from pathlib import Path

# Determine the base directory of the current file.
BASE_DIR: Path = Path(__file__).resolve().parent
# Construct the path to the 'data' directory relative to the base directory.
DATA_DIR: Path = BASE_DIR / "data"

# Load Transmembrane Domain (TMD) data from a JSON file.
with open(DATA_DIR / "tmd/tmd_list.json", "r") as f:
    TMD_DATA: dict = dict(json.load(f))

# Load Auto-inhibitory Peptide (AIP) data from a JSON file.
with open(DATA_DIR / "aip/aip_list.json", "r") as f:
    AIP_DATA: dict = dict(json.load(f))

# Load FRET-related Intracellular Domains (ICDs) data from a JSON file.
with open(DATA_DIR / "FRET/ICDs.json", "r") as f:
    FRET_ICDs: dict = dict(json.load(f))

# Load C-terminal TEV (Tobacco Etch Virus) protease chain data from a JSON file.
with open(DATA_DIR / "intracellular/CTEV_list.json", "r") as f:
    CTEV_DATA: dict = dict(json.load(f))

# Load N-terminal TEV (Tobacco Etch Virus) protease chain data from a JSON file.
with open(DATA_DIR / "intracellular/NTEV_list.json", "r") as f:
    NTEV_DATA: dict = dict(json.load(f))

# Load TEV Protease (TEVp) data from a JSON file.
with open(DATA_DIR / "intracellular/TEVp_list.json", "r") as f:
    TEVP_DATA: dict = dict(json.load(f))

# Load Protease Recognition Site (PRS) data from a JSON file.
with open(DATA_DIR / "prs/prs_list.json", "r") as f:
    PRS_DATA: dict = dict(json.load(f))

# Load signal sequences data from a JSON file.
with open(DATA_DIR / "signal_seqs/signal_sequences.json", "r") as f:
    SIGNAL_SEQS: dict = dict(json.load(f))

# Load tag sequences data from a JSON file.
with open(DATA_DIR / "tags/tag_sequences.json", "r") as f:
    TAG_SEQS: dict = dict(json.load(f))

# export all data for an easy overview
ALL_DATA = {
    "tmd": TMD_DATA,
    "aip": AIP_DATA,
    "fret": FRET_ICDs,
    "ctev": CTEV_DATA,
    "ntev": NTEV_DATA,
    "tev": TEVP_DATA,
    "prs": PRS_DATA,
    "signal": SIGNAL_SEQS,
    "tag": TAG_SEQS
}