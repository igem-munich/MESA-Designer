import sqlite3
import requests as req
from tqdm import tqdm
from pathlib import Path
from util.database_interaction import create_database, read_csv, create_table_from_header, insert_data, create_connection
import zipfile
import tarfile
import os

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

# create dictionaries for databases
Path.mkdir(Path("./data"))
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

# download the pdb files from abYbank's antibody database
print("Downloading abYbank Antibody DB (abdb)...")
if not download_file("http://www.abybank.org/abdb/Data/LH_Combined_Kabat.tar.bz2", "./files/abdb_pdbs.tar.bz2"):
    print("Failed to download abdb structures! Exiting!")
    exit(1)

# extract the tar archive and move files to the correct directory
with tarfile.open("./files/abdb_pdbs.tar.bz2", "r") as f:
    Path.mkdir(Path("./files/abdb_structures"), exist_ok=True)
    f.extractall("./files/abdb_structures")

# rename pdb files to be in line with sabdab naming (only chothia numbering)
for path in Path("./files/abdb_structures/chothia").glob("*.cho"):
    new_path = path.with_name(path.name[3:-3] + "pdb")
    os.rename(path, new_path)

print("Successfully downloaded abdb pdb files!")

# download pdb files from skempi database
print("Downloading SKEMPI structures...")
if not download_file("https://life.bsc.es/pid/skempi2/database/download/SKEMPI2_PDBs.tgz", "./files/skempi_pdbs.tgz"):
    print("Failed to download skempi structures! Exiting!")
    exit(1)

# extract tar archive and move files to correct directory
with tarfile.open("./files/skempi_pdbs.tgz", "r") as f:
    Path.mkdir(Path("./files/skempi_structures"), exist_ok=True)
    f.extractall("./files/skempi_structures")

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

# download sabdab pdb files and extract to correct location
print("Downloading sabdab pdb files...")
if not download_file("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/archive/all/", "./files/sabdab_structures.zip"):
    print("Failed to download sabdab structures! Exiting!")
    exit(1)

with zipfile.ZipFile("./files/sabdab_structures.zip", "r") as f:
    Path.mkdir(Path("./files/sabdab_structures"), exist_ok=True)
    f.extractall("./files/")

print("Successfully downloaded sabdab pdb files!")

# TODO: (optional) clean-up and tracking of installed files