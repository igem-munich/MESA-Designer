from pathlib import Path
import json

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR.parent / "data"
FILES_DIR = BASE_DIR.parent / "files"
RESOURCES_DIR = BASE_DIR.parent / "resources"

with open(RESOURCES_DIR / "tmd/tmd_list.json", "r") as f:
    TMD_DATA: dict = dict(json.load(f))