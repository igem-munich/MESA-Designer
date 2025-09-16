from pathlib import Path
import json

BASE_DIR: Path = Path(__file__).resolve().parent
DATA_DIR: Path = BASE_DIR.parent / "data"
FILES_DIR: Path = BASE_DIR.parent / "files"
RESOURCES_DIR: Path = BASE_DIR.parent / "resources"

with open(RESOURCES_DIR / "tmd/tmd_list.json", "r") as f:
    TMD_DATA: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "aip/aip_list.json", "r") as f:
    AIP_DATA: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "FRET/ICDs.json", "r") as f:
    FRET_ICDs: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "intracellular/CTEV_list.json", "r") as f:
    CTEV_DATA: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "intracellular/NTEV_list.json", "r") as f:
    NTEV_DATA: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "intracellular/TEVp_list.json", "r") as f:
    TEVP_DATA: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "prs/prs_list.json", "r") as f:
    PRS_DATA: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "colors/chain_colors.json", "r") as f:
    CHAIN_COLORS: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "signal_seqs/signal_sequences.json", "r") as f:
    SIGNAL_SEQS: dict[str, list[str]] = dict(json.load(f))

with open(RESOURCES_DIR / "tags/tag_sequences.json", "r") as f:
    TAG_SEQS: dict[str, list[str]] = dict(json.load(f))