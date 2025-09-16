from util import DATA_DIR, FILES_DIR
from util.database_interaction import *
from pathlib import Path
import datetime
import pandas as pd

conn: sqlite3.Connection = create_connection(str(DATA_DIR / "sabdab_summary_all.sqlite"))
sabdab_df: pd.DataFrame = get_dataframe(conn, "main")
skempi_df: pd.DataFrame = get_dataframe(conn, "skempi")

file_priority_list: list[str] = ["imgt", "chothia", "raw", "skempi", "abdb"]
skempi_pdbs: set[str] = set([x[:5] for x in list(skempi_df["#Pdb"])])

def get_highest_priority_path(list_of_paths: list[Path], priority_list: list[str]) -> Path | None:
    priority_map: dict[str, int] = {keyword: i for i, keyword in enumerate(priority_list)}

    def get_file_priority(file_path: Path) -> int:
        for keyword, priority in priority_map.items():
            if keyword in file_path.name.lower():
                return priority
        return len(priority_list)

    # Sort the list of paths using the custom key
    sorted_paths: list[Path] = sorted(list_of_paths, key=get_file_priority)

    # Return the first element, which is the highest-priority path
    return sorted_paths[0] if sorted_paths else None

def search_antibodies(antigen: str, filter_structures: bool=True) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], dict[str, str | Path], datetime.timedelta]:
    time: datetime.datetime = datetime.datetime.now()
    sabdab_selection: pd.DataFrame = sabdab_df.loc[sabdab_df["antigen_name"].str.contains(antigen, case=False) |
                                     sabdab_df["compound"].str.contains(antigen, case=False)].sort_values("affinity")
    skempi_selection: dict[str, pd.DataFrame] = {}
    pdb_files: dict[str, str | Path | list[Path]] = {}

    for _, row in sabdab_selection.iterrows():
        pdb: str = row["pdb"]

        # retrieve pdb files
        pdb_files[pdb] = []

        if Path(str(FILES_DIR) + "/sabdab_structures/imgt/" + pdb + ".pdb").is_file():
            pdb_files[pdb] = [Path(str(FILES_DIR) + "/sabdab_structures/imgt/" + pdb + ".pdb")]
        if not pdb_files[pdb]:
            if Path(str(FILES_DIR) + "/sabdab_structures/chothia/" + pdb + ".pdb").is_file():
                pdb_files[pdb] = [Path(str(FILES_DIR) + "/sabdab_structures/chothia/" + pdb + ".pdb")]
        if not pdb_files[pdb]:
            if Path(str(FILES_DIR) + "/sabdab_structures/raw/" + pdb + ".pdb").is_file():
                pdb_files[pdb] = [Path(str(FILES_DIR) + "/sabdab_structures/raw/" + pdb + ".pdb")]
        if not pdb_files[pdb]:
            pdb_files[pdb] = get_pdbs(pdb, str(FILES_DIR) + "/skempi_structures")
        if not pdb_files[pdb]:
            pdb_files[pdb] = get_pdbs(pdb, str(FILES_DIR) + "/abdb_structures/chothia")

        if len(pdb_files[pdb]) == 0:
            pdb_files[pdb] = ""
            continue
        pdb_files[pdb] = pdb_files[pdb][0]

        if pdb not in skempi_pdbs:
            continue

        # retrieve skempi data
        sel: pd.DataFrame = skempi_df.loc[skempi_df["#Pdb"].str[:5] == pdb]
        skempi_selection[pdb] = sel if len(sel) > 0 else None

    #    if filter_structures:
#        for key in pdb_files:
#            pdb_files[key] = get_highest_priority_path(pdb_files[key], file_priority_list)
    return sabdab_selection, skempi_selection, pdb_files, datetime.datetime.now()-time

def search_antibodies_api(antigen: str) -> dict[str, list[dict[str, str]] | float]:
    sabdab_selection, skempi_selection, pdb_files, search_duration = search_antibodies(antigen)

    # convert sabdab to dict
    sabdab_data = sabdab_selection.to_dict(orient="records")

    return {"sabdab_data": sabdab_data, "search_duration": search_duration.total_seconds()}