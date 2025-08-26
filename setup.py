import sqlite3
import requests as req
from tqdm import tqdm
from pathlib import Path
from util.database_interaction import create_database, read_csv, create_table_from_header, insert_data, create_connection

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
    #    RequiredFile("sabdab_sqlite", "./data/sabdab_summary_all.sqlite", FileType.DATABASE),
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
    header, data = read_csv("./files/sabdab_summary_all.tsv", delimiter="\t")
    create_table_from_header(conn, header, "main")
    insert_data(conn, header, "main")
    conn.close()

