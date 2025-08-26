import sqlite3
import numpy as np
import pandas as pd
import pathlib
import csv

# create new sqlite database
def create_databse(filepath: str) -> sqlite3.Connection:
    return sqlite3.Connection(filepath)

# create new table and populate with data
def read_csv(filepath: str, delimiter: str=","):
    data = []
    with open(filepath, 'r', newline='', encoding='utf-8') as csvfile:
        # Use csv.reader with tab as the delimiter
        reader = csv.reader(csvfile, delimiter=delimiter)
        header = next(reader) # Get the header row
        for row in reader:
            data.append(row)
    return header, data

def create_table_from_header(conn: sqlite3.Connection, header: str, table_name: str):
    columns = [f"'{col}' TEXT" for col in header]
    column_str = ', '.join(columns)

    create_table = f"create table {table_name} ({column_str})"

    try:
        cursor = conn.cursor()
        cursor.execute(create_table)
        conn.commit()
        print(f"Successfully created table {table_name} with columns: {column_str}")
    except sqlite3.Error as e:
        print(e)

def insert_data(conn: sqlite3.Connection, data: list, table_name: str):
    placeholders = ", ".join(["?" for _ in range(len(data[0]))])
    sql_insert = f"insert into {table_name} values ({placeholders})"

    try:
        cursor = conn.cursor()
        cursor.executemany(sql_insert, data)
        conn.commit()
        print("Data inserted successfully!")
    except sqlite3.Error as e:
        print(e)

# interact with existing databases
def create_connection(db_file: str):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to database: {db_file}")
    except sqlite3.Error as e:
        print(e)
    return conn

def retrieve_columns(conn: sqlite3.Connection, table_name: str):
    cursor = conn.cursor()
    cursor.execute(f"select * from {table_name}")
    columns = [desc[0] for desc in cursor.description]
    cursor.close()
    return columns

def get_entries(conn: sqlite3.Connection, table_name: str, key: str, value: str):
    valid_columns = retrieve_columns(conn, table_name)

    if not key in valid_columns:
        raise Exception(f"Illegal key, available columns are: {valid_columns}")

    cursor = conn.cursor()
    cursor.execute(f"select * from {table_name} where {key}=?", (value, ))
    rows = cursor.fetchall()
    cursor.close()
    return rows

def get_dataframe(conn: sqlite3.Connection, table_name: str):
    cursor = conn.cursor()
    cursor.execute(f"select * from {table_name}")
    rows = cursor.fetchall()
    columns = retrieve_columns(conn, table_name)
    arr = np.array([list(entry) for entry in rows])
    df = pd.DataFrame(arr, columns=columns)
    df = df.map(lambda cell : None if cell=="None" else cell)
    return df

def get_pdbs(pdb_id: str, directory: str):
    files = [file for file in pathlib.Path(directory).rglob("*.pdb") if str(file).__contains__(pdb_id)]
    return files

def search_antigen(conn: sqlite3.Connection, antigen_name: str):
    cursor = conn.cursor()
    cursor.execute(f"select * from main where antigen_name like ?", (f"%{antigen_name.lower()}%", ))
    rows = cursor.fetchall()
    return [list(row) for row in rows]

def extract_pdb_from_skempi(skempi_entry: str):
    return skempi_entry.split("_")[0]

def extract_pdb_from_abdb(abdb_entry: str):
    return abdb_entry[3:].split("_")[0]