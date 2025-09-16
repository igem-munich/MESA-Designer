import sqlite3
import numpy as np
import pandas as pd
import pathlib
import csv

# create new sqlite database
def create_database(filepath: str) -> sqlite3.Connection:
    return sqlite3.connect(filepath)

# create new table and populate with data
def read_csv(filepath: str, delimiter: str=",") -> tuple[list[str], list[list[str]]]:
    data: list[list[str]] = []
    with open(filepath, 'r', newline='', encoding='utf-8') as csvfile:
        # Use csv.reader with tab as the delimiter
        reader = csv.reader(csvfile, delimiter=delimiter)
        header: list[str] = next(reader) # Get the header row
        for row in reader:
            data.append(row)
    return header, data

def create_table_from_header(conn: sqlite3.Connection, header: list[str], table_name: str) -> None:
    columns: list[str] = [f"'{col}' TEXT" for col in header]
    column_str: str = ', '.join(columns)

    create_table: str = f"create table {table_name} ({column_str})"

    try:
        cursor: sqlite3.Cursor = conn.cursor()
        cursor.execute(create_table)
        conn.commit()
        print(f"Successfully created table {table_name} with columns: {column_str}")
    except sqlite3.Error as e:
        print(e)

def insert_data(conn: sqlite3.Connection, data: list, table_name: str) -> None:
    placeholders: str = ", ".join(["?" for _ in range(len(data[0]))])
    sql_insert: str = f"insert into {table_name} values ({placeholders})"

    try:
        cursor: sqlite3.Cursor = conn.cursor()
        cursor.executemany(sql_insert, data)
        conn.commit()
        print("Data inserted successfully!")
    except sqlite3.Error as e:
        print(e)

# interact with existing databases
def create_connection(db_file: str) -> sqlite3.Connection | None:
    conn: sqlite3.Connection | None = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to database: {db_file}")
    except sqlite3.Error as e:
        print(e)
    return conn

def retrieve_columns(conn: sqlite3.Connection, table_name: str) -> list[str]:
    cursor: sqlite3.Cursor = conn.cursor()
    cursor.execute(f"select * from {table_name}")
    columns: list[str] = [desc[0] for desc in cursor.description]
    cursor.close()
    return columns

def get_entries(conn: sqlite3.Connection, table_name: str, key: str, value: str) -> list[list[str]]:
    valid_columns: list[str] = retrieve_columns(conn, table_name)

    if not key in valid_columns:
        raise Exception(f"Illegal key, available columns are: {valid_columns}")

    cursor: sqlite3.Cursor = conn.cursor()
    cursor.execute(f"select * from {table_name} where {key}=?", (value, ))
    rows: list[list[str]] = cursor.fetchall()
    cursor.close()
    return rows

def get_dataframe(conn: sqlite3.Connection, table_name: str) -> pd.DataFrame:
    cursor: sqlite3.Cursor = conn.cursor()
    cursor.execute(f"select * from {table_name}")
    rows: list[list[str]] = cursor.fetchall()
    columns: list[str] = retrieve_columns(conn, table_name)
    arr: np.ndarray = np.array([list(entry) for entry in rows])
    df: pd.DataFrame = pd.DataFrame(arr, columns=columns)
    df = df.map(lambda cell : None if cell=="None" else cell)
    return df

def get_pdbs(pdb_id: str, directory: str) -> list[pathlib.Path]:
    files = [file for file in pathlib.Path(directory).rglob(f"{pdb_id}*.pdb")]
    return files

def search_antigen(conn: sqlite3.Connection, antigen_name: str) -> list[list[str]]:
    cursor: sqlite3.Cursor = conn.cursor()
    cursor.execute(f"select * from main where antigen_name like ?", (f"%{antigen_name.lower()}%", ))
    rows: list[list[str]] = cursor.fetchall()
    return [list(row) for row in rows]

def extract_pdb_from_skempi(skempi_entry: str) -> str:
    return skempi_entry.split("_")[0]

def extract_pdb_from_abdb(abdb_entry: str) -> str:
    return abdb_entry[3:].split("_")[0]