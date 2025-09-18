import sqlite3
import numpy as np
import pandas as pd
import pathlib
import csv


# create new sqlite database
def create_database(filepath: str) -> sqlite3.Connection:
    """
    Creates a new SQLite database file at the specified path and returns a connection object.
    :param filepath: The path where the SQLite database file will be created or opened.
    :return: A connection object to the SQLite database.
    """
    return sqlite3.connect(filepath)


# create new table and populate with data
def read_csv(filepath: str, delimiter: str=",") -> tuple[list[str], list[list[str]]]:
    """
    Reads data from a CSV file, including the header and all rows.
    :param filepath: The path to the CSV file.
    :param delimiter: The delimiter used in the CSV file (default is comma).
    :return: A tuple containing two elements: a list of strings representing the header, and a list of lists of strings representing the data rows.
    """
    data: list[list[str]] = []
    with open(filepath, 'r', newline='', encoding='utf-8') as csvfile:
        # Use csv.reader with tab as the delimiter
        reader = csv.reader(csvfile, delimiter=delimiter)
        header: list[str] = next(reader) # Get the header row
        for row in reader:
            data.append(row)
    return header, data


def create_table_from_header(conn: sqlite3.Connection, header: list[str], table_name: str) -> None:
    """
    Creates a new table in the specified SQLite database using a list of column headers.
    All columns are created with TEXT data type.
    :param conn: The SQLite database connection object.
    :param header: A list of strings, where each string is a column name for the new table.
    :param table_name: The name of the table to be created.
    :return: None
    """
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
    """
    Inserts multiple rows of data into a specified table in the SQLite database.
    :param conn: The SQLite database connection object.
    :param data: A list of lists, where each inner list represents a row of data to be inserted.
    :param table_name: The name of the table into which data will be inserted.
    :return: None
    """
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
    """
    Establishes a connection to an existing SQLite database file.
    :param db_file: The path to the SQLite database file.
    :return: An sqlite3.Connection object if the connection is successful, otherwise None.
    """
    conn: sqlite3.Connection | None = None
    try:
        conn = sqlite3.connect(db_file)
        print(f"Connected to database: {db_file}")
    except sqlite3.Error as e:
        print(e)
    return conn


def retrieve_columns(conn: sqlite3.Connection, table_name: str) -> list[str]:
    """
    Retrieves the column names of a specified table from the SQLite database.
    :param conn: The SQLite database connection object.
    :param table_name: The name of the table whose columns are to be retrieved.
    :return: A list of strings, where each string is a column name.
    """
    cursor: sqlite3.Cursor = conn.cursor()
    cursor.execute(f"select * from {table_name}")
    columns: list[str] = [desc[0] for desc in cursor.description]
    cursor.close()
    return columns


def get_entries(conn: sqlite3.Connection, table_name: str, key: str, value: str) -> list[list[str]]:
    """
    Retrieves entries from a specified table where a given column (key) matches a specific value.
    :param conn: The SQLite database connection object.
    :param table_name: The name of the table to query.
    :param key: The name of the column to use for filtering.
    :param value: The value to match against the specified column.
    :return: A list of lists of strings, where each inner list represents a matching row.
    :raises Exception: If the provided key is not a valid column in the table.
    """
    valid_columns: list[str] = retrieve_columns(conn, table_name)

    if not key in valid_columns:
        raise Exception(f"Illegal key, available columns are: {valid_columns}")

    cursor: sqlite3.Cursor = conn.cursor()
    cursor.execute(f"select * from {table_name} where {key}=?", (value, ))
    rows: list[list[str]] = cursor.fetchall()
    cursor.close()
    return rows


def get_dataframe(conn: sqlite3.Connection, table_name: str) -> pd.DataFrame:
    """
    Retrieves all data from a specified table and returns it as a pandas DataFrame.
    String "None" values in the database are converted to Python's None.
    :param conn: The SQLite database connection object.
    :param table_name: The name of the table to retrieve data from.
    :return: A pandas DataFrame containing all data from the specified table.
    """
    cursor: sqlite3.Cursor = conn.cursor()
    cursor.execute(f"select * from {table_name}")
    rows: list[list[str]] = cursor.fetchall()
    columns: list[str] = retrieve_columns(conn, table_name)
    arr: np.ndarray = np.array([list(entry) for entry in rows])
    df: pd.DataFrame = pd.DataFrame(arr, columns=columns)
    df = df.map(lambda cell : None if cell=="None" else cell)
    return df


def get_pdbs(pdb_id: str, directory: str) -> list[pathlib.Path]:
    """
    Searches for PDB files within a specified directory and its subdirectories that match a given PDB ID.
    :param pdb_id: The PDB ID to search for (e.g., "1abc").
    :param directory: The directory path where the search for PDB files will start.
    :return: A list of pathlib.Path objects for all matching PDB files.
    """
    files = [file for file in pathlib.Path(directory).rglob(f"{pdb_id}*.pdb")]
    return files


def search_antigen(conn: sqlite3.Connection, antigen_name: str) -> list[list[str]]:
    """
    Searches for entries in the 'main' table where the 'antigen_name' column contains the given antigen name (case-insensitive).
    :param conn: The SQLite database connection object.
    :param antigen_name: The name or partial name of the antigen to search for.
    :return: A list of lists of strings, where each inner list represents a matching row.
    """
    cursor: sqlite3.Cursor = conn.cursor()
    cursor.execute(f"select * from main where antigen_name like ?", (f"%{antigen_name.lower()}%", ))
    rows: list[list[str]] = cursor.fetchall()
    return [list(row) for row in rows]


def extract_pdb_from_skempi(skempi_entry: str) -> str:
    """
    Extracts the PDB ID from a SKEMPI entry string.
    It assumes the PDB ID is the part before the first underscore.
    :param skempi_entry: The SKEMPI entry string (e.g., "1ABC_A_B").
    :return: The extracted PDB ID (e.g., "1ABC").
    """
    return skempi_entry.split("_")[0]


def extract_pdb_from_abdb(abdb_entry: str) -> str:
    """
    Extracts the PDB ID from an ABDB entry string.
    It assumes the PDB ID starts from the 4th character and is the part before the first underscore.
    :param abdb_entry: The ABDB entry string (e.g., "ab_1ABC_chain").
    :return: The extracted PDB ID (e.g., "1ABC").
    """
    return abdb_entry[3:].split("_")[0]