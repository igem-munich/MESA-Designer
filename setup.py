import sqlite3
import requests as req
from tqdm import tqdm
from pathlib import Path
from util.database_interaction import create_database, read_csv, create_table_from_header, insert_data, create_connection
import zipfile
import tarfile
import os
import shutil

# Defines a function to download a file from a given URL to a specified file path
def download_file(url: str, file_path: str) -> bool:
    try:
        # Create parent directories for the file if they don't already exist
        Path.mkdir(Path(file_path).parent, parents=True, exist_ok=True)

        # Send a GET request to the URL, enabling streaming to handle large files
        res = req.get(url, stream=True)
        # Raise an HTTPError for bad responses (4xx or 5xx)
        res.raise_for_status()

        # Get the total size of the file from the 'content-length' header for the progress bar. Defaults to 0 for missing header
        size = int(res.headers.get("content-length", 0))

        # Initialize a tqdm progress bar to display download progress
        # It shows total size, unit (Bytes), unit scaling (to KB/MB), and a description
        with tqdm(total=size, unit="B", unit_scale=True, unit_divisor=1024, desc=file_path.split("/")[-1]) as progress:
            # Open the specified file in binary write mode
            with open(file_path, "wb") as f:
                # Iterate over the content in chunks and write to the file
                for chunk in res.iter_content(chunk_size=8192):
                    f.write(chunk)
                    # Update the progress bar with the size of the written chunk
                    progress.update(len(chunk))

        # Return True if the download was successful
        return True

    except Exception as e:
        # Catch any exceptions during download and print the error message.
        print(e)
        # Return False if an error occurred during download.
        return False

# check if download is complete
if Path("./files/done.txt").exists():
    print("Files already downloaded!\nIf you wish to re-download, please remove '/files/done.txt'")
    exit()

# online mode
offline_mode = False # input("Would you like to use this offline? Note: requires 55GB download

# create dictionaries for databases
Path.mkdir(Path("./data"))

if offline_mode:
    Path.mkdir(Path("./files"))

# Section for downloading the SABDAB database and setting up its SQLite database
print("Downloading sabdab database and creating sqlite database...")
print("Downloading summary file...")
# Attempt to download the SABDAB summary TSV file
if not download_file("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/", "./files/sabdab_summary_all.tsv"):
    print("Failed to download sabdab summary! Exiting!")
    exit(1)

# if successful, create a new sqlite3 database and necessary parent directories. then create a "main" table and insert all data from the csv file into the table
print("Creating sqlite3 database...")
Path.mkdir(Path("./data"), parents=True, exist_ok=True)
conn: sqlite3.Connection = create_database("./data/sabdab_summary_all.sqlite")
header, data = read_csv("./files/sabdab_summary_all.tsv", "\t")
create_table_from_header(conn, header, "main")
insert_data(conn, data, "main")
conn.close()
print("Successfully created sabdab sqlite3 database")

if offline_mode:
    # download pdb files from skempi database
    print("Downloading SKEMPI structures...")
    if not download_file("https://life.bsc.es/pid/skempi2/database/download/SKEMPI2_PDBs.tgz", "./files/skempi_pdbs.tgz"):
        print("Failed to download skempi structures! Exiting!")
        exit(1)

    # extract tar archive and move files to correct directory
    print("Extracting archive...")
    with tarfile.open("./files/skempi_pdbs.tgz", "r") as f:
        Path.mkdir(Path("./files/skempi_structures_tmp"), exist_ok=True)
        f.extractall("./files/skempi_structures_tmp")

    # move pdb files to correct location
    skempi_dir: Path = Path("./files/skempi_structures_tmp/PDBs")
    shutil.move(skempi_dir, skempi_dir.parent.parent)
    shutil.rmtree(skempi_dir.parent)
    os.rename(Path("./files/PDBs"), Path("./files/skempi_structures"))

    # remove all files which start with ._
    for path in Path("./files/skempi_structures").glob("._*"):
        os.remove(path)

    # rename pdb files to be in line with sabdab naming (entire database)
    for path in Path("./files/skempi_structures").glob("*"):
        new_path = path.with_name(path.name.lower())
        os.rename(path, new_path)

    print("Successfully downloaded skempi pdb files!")

# download skempi protein-protein affinity data and insert it into sqlite database
print("Downloading skempi affinity data...")
if not download_file("https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv", "./files/skempi_v2.csv"):
    print("Failed to download skempi v2 affinity data! Exiting!")
    exit(1)

# create a new table and insert the data
print("Adding data to sqlite database...")
skempi_header, skempi_data = read_csv("./files/skempi_v2.csv", ";")
conn = create_connection("./data/sabdab_summary_all.sqlite")
create_table_from_header(conn, skempi_header, "skempi")
insert_data(conn, skempi_data, "skempi")
print("Successfully inserted skempi data into database!")

if offline_mode:
    # download the pdb files from abYbank's antibody database
    print("Downloading abYbank Antibody DB (abdb)...")
    if not download_file("http://www.abybank.org/abdb/snapshots/abdb_20240706.zip", "./files/abdb_pdbs.zip"):
        print("Failed to download abdb structures! Exiting!")
        exit(1)

    # extract the tar archive and move files to the correct directory
    print("Extracting archive...")
    with zipfile.ZipFile("./files/abdb_pdbs.zip", "r") as f:
        Path.mkdir(Path("./files/abdb_structures_tmp"), exist_ok=True)
        f.extractall("./files/abdb_structures_tmp")

    # move content of downloaded snapshot to parent directory
    # get folder containing the files
    # needs to be updated when abdb is done rebuilding their database
    abdb_dir: Path = Path("./files/abdb_structures_tmp/abdb_newdata_20240706")
    shutil.move(abdb_dir, abdb_dir.parent.parent)
    shutil.rmtree(abdb_dir.parent)
    os.rename(Path("./files/abdb_newdata_20240706"), Path("./files/abdb_structures"))

    # move chothia files to separate directory
    Path.mkdir(Path("./files/abdb_structures/chothia"))
    for path in Path("./files/abdb_structures").glob("*.cho"):
        shutil.move(path, path.parent / "chothia")

    # rename pdb files to be in line with sabdab naming (only chothia numbering)
    for path in Path("./files/abdb_structures/chothia").glob("*.cho"):
        new_path = path.with_name(path.name[3:-3] + "pdb")
        os.rename(path, new_path)

    print("Successfully downloaded abdb pdb files!")

if offline_mode:
    # download sabdab pdb files and extract to correct location
    print("Downloading sabdab pdb files...")
    if not download_file("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/archive/all/", "./files/sabdab_structures.zip"):
        print("Failed to download sabdab structures! Exiting!")
        exit(1)

    # extract zip archive and move to correct location
    print("Extracting archive...")
    with zipfile.ZipFile("./files/sabdab_structures.zip", "r") as f:
        Path.mkdir(Path("./files/sabdab_structures_tmp"), exist_ok=True)
        f.extractall("./files/sabdab_structures_tmp")

    sabdab_dir: Path = Path("./files/sabdab_structures_tmp/all_structures")
    shutil.move(sabdab_dir, sabdab_dir.parent.parent)
    shutil.rmtree(sabdab_dir.parent)
    os.rename(Path("./files/all_structures"), Path("./files/sabdab_structures"))

    print("Successfully downloaded sabdab pdb files!")

print("Successfully setup databases!")

# track completed download
open("./files/done.txt", "a").close()

# TODO: (optional) clean-up and tracking of installed files