from util import DATA_DIR, FILES_DIR
from util.database_interaction import *
from pathlib import Path
import datetime

conn = create_connection(str(DATA_DIR / "sabdab_summary_all.sqlite"))
sabdab_df = get_dataframe(conn, "main")
skempi_df = get_dataframe(conn, "skempi")

file_priority_list = ["imgt", "chothia", "raw", "skempi", "abdb"]
skempi_pdbs = set([x[:5] for x in list(skempi_df["#Pdb"])])

def get_highest_priority_path(list_of_paths, priority_list):
    priority_map = {keyword: i for i, keyword in enumerate(priority_list)}

    def get_file_priority(file_path: Path):
        for keyword, priority in priority_map.items():
            if keyword in file_path.name.lower():
                return priority
        return len(priority_list)

    # Sort the list of paths using the custom key
    sorted_paths = sorted(list_of_paths, key=get_file_priority)

    # Return the first element, which is the highest-priority path
    return sorted_paths[0] if sorted_paths else None

def search_antibodies(antigen: str, filter_structures: bool=True):
    time = datetime.datetime.now()
    sabdab_selection = sabdab_df.loc[sabdab_df["antigen_name"].str.contains(antigen, case=False) |
                                     sabdab_df["compound"].str.contains(antigen, case=False)].sort_values("affinity")
    skempi_selection = {}
    pdb_files = {}

    for _, row in sabdab_selection.iterrows():
        pdb = row["pdb"]
        
        # retrieve pdb files
        pdb_files[pdb] = []
        pdb_files[pdb] = get_pdbs(pdb, str(FILES_DIR)+"/abdb_structures/chothia")
        if pdb_files[pdb] == []:
            pdb_files[pdb] = get_pdbs(pdb, str(FILES_DIR)+"/skempi_structures")
        if pdb_files[pdb] == []:
            pdb_files[pdb] = get_pdbs(pdb, str(FILES_DIR)+"/sabdab_structures")
        
        if(len(pdb_files[pdb]) == 0):
            pdb_files[pdb] = ""
            continue
        pdb_files[pdb] = pdb_files[pdb][0]

        if pdb not in skempi_pdbs:
            continue
        
        # retrieve skempi data
        sel = skempi_df.loc[skempi_df["#Pdb"].str[:5]==pdb]
        skempi_selection[pdb] = sel if len(sel) > 0 else None

        
#    if filter_structures:
#        for key in pdb_files:
#            pdb_files[key] = get_highest_priority_path(pdb_files[key], file_priority_list)
    return sabdab_selection, skempi_selection, pdb_files, datetime.datetime.now()-time