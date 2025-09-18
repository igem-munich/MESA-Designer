from util import DATA_DIR, FILES_DIR
from util.database_interaction import *
from pathlib import Path
import datetime
import pandas as pd

conn: sqlite3.Connection = create_connection(str(DATA_DIR / "sabdab_summary_all.sqlite")) # Establishes a connection to an SQLite database file named "sabdab_summary_all.sqlite" located in DATA_DIR.
sabdab_df: pd.DataFrame = get_dataframe(conn, "main") # Loads data from the "main" table of the connected SQLite database into a pandas DataFrame named sabdab_df.
skempi_df: pd.DataFrame = get_dataframe(conn, "skempi") # Loads data from the "skempi" table of the connected SQLite database into a pandas DataFrame named skempi_df.

file_priority_list: list[str] = ["imgt", "chothia", "raw", "skempi", "abdb"] # Defines a list of keywords representing the priority order for different types of PDB files.
skempi_pdbs: set[str] = set([x[:5] for x in list(skempi_df["#Pdb"])]) # Strips extra annotation of skempi pdb files to be compatible with regular pdb ids.

def get_highest_priority_path(list_of_paths: list[Path], priority_list: list[str]) -> Path | None: # Defines a function to get the highest priority file path from a list, based on a given priority list of keywords.
    priority_map: dict[str, int] = {keyword: i for i, keyword in enumerate(priority_list)} # Creates a dictionary mapping keywords from the priority_list to their integer priorities (lower number means higher priority).

    def get_file_priority(file_path: Path) -> int: # Defines a nested helper function to determine the priority of a given file path.
        for keyword, priority in priority_map.items(): # Iterates through the keyword-priority pairs in the priority_map.
            if keyword in file_path.name.lower(): # Checks if the keyword is present in the lowercase version of the file path. #TODO: might be a bug (checks name instead of path)
                return priority # Returns the priority if a matching keyword is found.
        return len(priority_list) # Returns a low priority (length of the list) if no keyword matches.

    sorted_paths: list[Path] = sorted(list_of_paths, key=get_file_priority) # Sorts the input list of paths based on their assigned priorities using the get_file_priority function.

    return sorted_paths[0] if sorted_paths else None # Returns the first path in the sorted list (highest priority), or None if the list is empty.

def search_antibodies(antigen: str, filter_structures: bool=True) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], dict[str, str | Path], datetime.timedelta]: # Defines a function to search for antibodies based on an antigen, with an option to filter structures.
    time: datetime.datetime = datetime.datetime.now() # Records the current time to measure the search duration.
    sabdab_selection: pd.DataFrame = sabdab_df.loc[sabdab_df["antigen_name"].str.contains(antigen, case=False) | # Selects rows from sabdab_df where the 'antigen_name' column (case-insensitive)
                                     sabdab_df["compound"].str.contains(antigen, case=False)].sort_values("affinity") # or 'compound' column contains the given antigen, then sorts by 'affinity'.
    skempi_selection: dict[str, pd.DataFrame] = {} # Initializes an empty dictionary to store SKEMPI data related to selected PDBs.
    pdb_files: dict[str, str | Path | list[Path]] = {} # Initializes an empty dictionary to store file paths for PDBs.

    for _, row in sabdab_selection.iterrows(): # Iterates through each row in the sabdab_selection DataFrame.
        pdb: str = row["pdb"] # Extracts the PDB ID from the current row.

        pdb_files[pdb] = [] # Initializes an empty list for the current PDB ID to store potential file paths.

        if Path(str(FILES_DIR) + "/sabdab_structures/imgt/" + pdb + ".pdb").is_file(): # Checks if an IMGT-formatted PDB file exists for the current PDB ID.
            pdb_files[pdb] = [Path(str(FILES_DIR) + "/sabdab_structures/imgt/" + pdb + ".pdb")] # If found, adds the IMGT PDB file path to the list.
        if not pdb_files[pdb]: # If no IMGT PDB file was found.
            if Path(str(FILES_DIR) + "/sabdab_structures/chothia/" + pdb + ".pdb").is_file(): # Checks if a Chothia-formatted PDB file exists.
                pdb_files[pdb] = [Path(str(FILES_DIR) + "/sabdab_structures/chothia/" + pdb + ".pdb")] # If found, adds the Chothia PDB file path.
        if not pdb_files[pdb]: # If neither IMGT nor Chothia PDB file was found.
            if Path(str(FILES_DIR) + "/sabdab_structures/raw/" + pdb + ".pdb").is_file(): # Checks if a raw PDB file exists.
                pdb_files[pdb] = [Path(str(FILES_DIR) + "/sabdab_structures/raw/" + pdb + ".pdb")] # If found, adds the raw PDB file path.
        if not pdb_files[pdb]: # If no PDB file was found in SABDAB structures.
            pdb_files[pdb] = get_pdbs(pdb, str(FILES_DIR) + "/skempi_structures") # Tries to retrieve PDB files from SKEMPI structures.
        if not pdb_files[pdb]: # If still no PDB file was found.
            pdb_files[pdb] = get_pdbs(pdb, str(FILES_DIR) + "/abdb_structures/chothia") # Tries to retrieve PDB files from ABDB Chothia structures.

        if len(pdb_files[pdb]) == 0: # If no PDB files were found after checking all sources.
            pdb_files[pdb] = "" # Sets the entry for the current PDB to an empty string.
            continue # Skips to the next iteration of the loop.
        pdb_files[pdb] = pdb_files[pdb][0] # If multiple PDB files were found, takes the first one (implicitly, the highest priority from 'get_pdbs').

        if pdb not in skempi_pdbs: # Checks if the current PDB ID is not present in the set of SKEMPI PDBs.
            continue # If not, skips to the next iteration as there won't be SKEMPI data.

        sel: pd.DataFrame = skempi_df.loc[skempi_df["#Pdb"].str[:5] == pdb] # Selects rows from skempi_df where the PDB ID matches the current one.
        skempi_selection[pdb] = sel if len(sel) > 0 else None # Stores the selected SKEMPI data for the PDB, or None if no data is found.

    #    if filter_structures:
#        for key in pdb_files:
#            pdb_files[key] = get_highest_priority_path(pdb_files[key], file_priority_list)
    return sabdab_selection, skempi_selection, pdb_files, datetime.datetime.now()-time # Returns the SABDAB selection, SKEMPI selection, PDB file paths, and the duration of the search.

def search_antibodies_api(antigen: str) -> dict[str, list[dict[str, str]] | float]: # Defines a function for searching antibodies, intended for API use, returning results as a dictionary.
    sabdab_selection, skempi_selection, pdb_files, search_duration = search_antibodies(antigen) # Calls the main search_antibodies function to get the results.

    # convert sabdab to dict
    sabdab_data = sabdab_selection.to_dict(orient="records") # Converts the sabdab_selection DataFrame into a list of dictionaries, where each dictionary represents a row.

    return {"sabdab_data": sabdab_data, "search_duration": search_duration.total_seconds()} # Returns a dictionary containing the SABDAB data and the search duration in seconds.