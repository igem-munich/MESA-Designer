import sqlite3
import requests as req
from tqdm import tqdm
from pathlib import Path
from util.database_interaction import create_database, read_csv, create_table_from_header, insert_data, create_connection
import zipfile

def download_file(url: str, file_path: str) -> bool:
    try:
        Path.mkdir(Path(file_path).parent, parents=True, exist_ok=True)

        res = req.get(url, stream=True)
        res.raise_for_status()

        size = int(res.headers.get("content-length", 0))

        with tqdm(total=size, unit="B", unit_scale=True, unit_divisor=1024, desc=file_path.split("/")[-1]) as progress:
            with open(file_path, "wb") as f:
                for chunk in res.iter_content(chunk_size=8192):
                    f.write(chunk)
                    progress.update(len(chunk))

        return True

    except Exception as e:
        print(e)
        return False

if  __name__ == "__main__":
    # constants
    #required_files: list[RequiredFile] = [
    #    RequiredFile("abdb_pdbs", "./abdb_pdbs", FileType.DIRECTORY),
    #    RequiredFile("sabdab_structures", "./sabdab_structures", FileType.DIRECTORY),
    #    RequiredFile("skempi2_pdbs", "./SKEMPI2_PDBS", FileType.DIRECTORY)
    #]
    Path.mkdir(Path("./data"))
    Path.mkdir(Path("./files"))

    # sabdab_sqlite
    print("Downloading sabdab database and creating sqlite database...")
    print("Downloading summary file...")
    if not download_file("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/", "./files/sabdab_summary_all.tsv"):
        print("Failed to download sabdab summary! Exiting!")
        exit(1)

    print("Creating sqlite3 database...")
    Path.mkdir(Path("./data"), parents=True, exist_ok=True)
    conn: sqlite3.Connection = create_database("./data/sabdab_summary_all.sqlite")
    header, data = read_csv("./files/sabdab_summary_all.tsv", "\t")
    create_table_from_header(conn, header, "main")
    insert_data(conn, data, "main")
    conn.close()
    print("Successfully created sabdab sqlite3 database")

    print("Downloading sabdab pdb files...")
    if not download_file("https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/archive/all/", "./files/sabdab_structures.zip"):
        print("Failed to download sabdab structures! Exiting!")
        exit(1)

    with zipfile.ZipFile("./files/sabdab_structures.zip", "r") as f:
        Path.mkdir(Path("./files/sabdab_structures"), exist_ok=True)
        f.extractall("./files/")
    print("Successfully downloaded sabdab pdb files!")

    print("Downloading abYbank Antibody DB (abdb)...")
    if not download_file("http://www.abybank.org/abdb/Data/LH_Combined_Kabat.tar.bz2", "./files/abdb_pdbs.tar.bz2"):
        print("Failed to download abdb structures! Exiting!")
        exit(1)

    with zipfile.ZipFile("./files/abdb_pdbs.tar.bz2", "r") as f:
        Path.mkdir(Path("./files/abdb_structures"), exist_ok=True)
        f.extractall("./files/abdb_structures")

    # TODO: rename pdbs according to renaming_to_pure_pdb.ipynb
    print("Successfully downloaded abdb pdb files!")

    print("Downloading SKEMPI structures...")
    if not download_file("https://life.bsc.es/pid/skempi2/database/download/SKEMPI2_PDBs.tgz", "./files/skempi_pdbs.tgz"):
        Path.mkdir(Path("./files/skempi_structures"), exist_ok=True)
        f.extractall("./files/skempi_structures")

    # TODO: rename pdbs according to renaming_to_pure_pdb.ipynb
    print("Successfully downloaded skempi pdb files!")

    # TODO: (optional) clean-up and tracking of installed files

    print("Downloading skempi affinity data...")
    if not download_file("https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv", "./files/skempi_v2.csv"):
        print("Failed to download skempi v2 affinity data! Exiting!")
        exit(1)

    print("Adding data to sqlite database...")
    skempi_header, skempi_data = read_csv("./files/skempi_v2.csv", ";")
    conn = create_connection("./data/sabdab_summary_all.sqlite")
    create_table_from_header(conn, skempi_header, "skempi")
    insert_data(conn, skempi_data, "skempi")
    print("Successfully inserted skempi data into database!")