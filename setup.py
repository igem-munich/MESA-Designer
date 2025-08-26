import requests as req
from enum import Enum, auto
from os import path
from typing import Final
from pathlib import Path

# create enum for different file types
class FileType(Enum):
    DIRECTORY = auto()
    FILE = auto()
    DATABASE = auto()

class RequiredFile:
    def __init__(self, name: str, file_path: str, filetype: FileType, url: str = None):
        self.file_path = file_path
        self.name = name
        self.filetype = filetype
        self.url = url

    def check_existence(self) -> bool:
        match self.filetype:
            case FileType.FILE:
                return path.exists(self.file_path)
            case FileType.DIRECTORY:
                return path.exists(self.file_path)
            case FileType.DATABASE:
                return path.exists(self.file_path)

    def get_file(self) -> bool:
        if self.filetype == FileType.DATABASE:
            raise Exception("Database needs to be created using create_database_from_csv()")




    def create_database_from_csv(self, file_path: str) -> bool:

    def __str__(self) -> str:
        return self.name

# download file to path
def get_file(url: str, file_path: str) -> bool:
    try:
        with req.get(url, stream=True) as r:
            r.raise_for_status()
            Path.mkdir(Path(file_path), parents=True, exist_ok=True)

            with open(file_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

        return True

    except Exception as e:
        print(e.with_traceback(e.__traceback__))
        return False

if  __name__ == "__main__":
    # constants
    required_files: list[RequiredFile] = [
        RequiredFile("sabdab_sqlite", "./data/sabdab_summary_all.sqlite", FileType.DATABASE),
        RequiredFile("abdb_pdbs", "./abdb_pdbs", FileType.DIRECTORY),
        RequiredFile("sabdab_structures", "./sabdab_structures", FileType.DIRECTORY),
        RequiredFile("skempi2_pdbs", "./SKEMPI2_PDBS", FileType.DIRECTORY)
    ]

    missing_files = [file for file in required_files if not file.check_existence()]

    # update user
    if len(missing_files) == 0:
        print("No missing files found, assuming current version and quitting!")
        exit(0)

    print(f"Found {len(missing_files)} missing files:")
    for file in missing_files:
        print(f"- {file}")

    # automatic download
    if not input("Would you like to continue with automatic download? [y/n]: ").lower() == "y":
        exit(0)

    for file in missing_files:



