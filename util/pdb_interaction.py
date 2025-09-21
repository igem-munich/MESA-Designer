from Bio import SeqIO
from io import StringIO
import re
import requests


def extract_chains_from_pdb(file_path: str | None = None, file_content: str | None = None) -> dict[str, dict[str, str | int]]:
    """
    Extracts chain data (ID, chain ID, FASTA name, sequence) from a PDB file, either from a file path or direct content.
    :param file_path: The path to the PDB file.
    :param file_content: The content of the PDB file as a string.
    :return: A list of dictionaries, where each dictionary represents a chain with its 'id', 'chain_id', 'fasta_name', and 'sequence'.
    """
    try:
        chains_data: dict[str, dict[str, str]] = {}
        records = SeqIO.parse(StringIO(file_content) if file_content else file_path, "pdb-atom")
        for record in records:
            fasta_name: str = re.sub("[^a-zA-Z0-9]", "", record.id.replace(":", "_"))

            chains_data[record.annotations["chain"]] = {
                "id": record.id,
                "chain_id": record.annotations["chain"],
                "fasta_name": fasta_name,
                "sequence": str(record.seq),
                "start": record.annotations["start"],
                "end": record.annotations["end"]
            }

        return chains_data
    except Exception as e:
        raise(e)


def extract_chains_from_fasta(fasta: str) -> list[dict] | None:
    """
    Extracts chain data from a FASTA formatted string.
    :param fasta: A string containing FASTA formatted sequences.
    :return: A list of Bio.SeqRecord.Record objects, or None if an error occurs.
    """
    chains_data: list = []
    try:
        fasta_io: StringIO = StringIO(fasta)
        records = SeqIO.parse(fasta_io, "fasta")

        for record in records:
            chains_data.append(record)

        return chains_data

    except Exception as e:
        print(e)
        return None


def convert_chains_to_fasta_string(chains_data: list[dict[str, str]]) -> str:
    """
    Converts a list of chain dictionaries to a FASTA formatted string.
    :param chains_data: A list of dictionaries, where each dictionary represents a chain with 'fasta_name' and 'sequence' keys.
    :return: A string containing the FASTA formatted sequences of the chains.
    """
    fasta_string: StringIO = StringIO()
    for chain in chains_data:
        fasta_string.write(f">{chain['fasta_name']}\n{chain['sequence']}\n")
    return fasta_string.getvalue()


def extract_fasta_from_pdb(file_path: str) -> str:
    """
    Extracts FASTA formatted sequences from a PDB file.
    Original function, modified to use extract_chains_from_pdb for consistency.
    :param file_path: The path to the PDB file.
    :return: A string containing the FASTA formatted sequences extracted from the PDB file.
    """
    chains_data: list[dict[str, str]] = extract_chains_from_pdb(file_path)
    return convert_chains_to_fasta_string(chains_data)


def cut_protein(file_content: str) -> str:
    """
    Function for cutting a protein from its file content.
    :param file_content: The content of the protein file as a string.
    :return: Raises an exception as the feature is currently under construction.
    """
    raise Exception("Feature currently under construction!")


def get_pdb_from_rcsb(pdb_id: str) -> str | None:
    """
    Fetches a PDB file in plain text format from the RCSB PDB database.
    :param pdb_id: The PDB ID of the structure to fetch.
    :return: The content of the PDB file as a string if successful, otherwise None.
    """
    url: str = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        res: requests.Response = requests.get(url)
        res.raise_for_status()
        return res.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching PDB file for {pdb_id}: {e}")
        return None


def get_fasta_from_rcsb(pdb_id: str) -> str | None:
    """
    Fetches FASTA formatted sequences for a given PDB ID from the RCSB PDB database.
    :param pdb_id: The PDB ID of the structure to fetch FASTA for.
    :return: The content of the FASTA file as a string if successful, otherwise None.
    """
    url: str = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    try:
        res: requests.Response = requests.get(url)
        res.raise_for_status()
        return res.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching FASTA file for {pdb_id}: {e}")
        return None


def generate_chain_selection(pdb_content: str, selection: dict[str, tuple[int, int]]) -> dict[str, str] | None:
    """
    Generates a FASTA string from selected chains and residues within PDB content.
    :param pdb_content: The full pdb file as a string.
    :param selection: A dictionary where keys are chain IDs and values are tuples representing the (start_residue_index, end_residue_index) (0-indexed, end exclusive).
    :return: A dictionary where keys are chain IDs and values are the selected sequence strings, or None if input is invalid.
    """
    if not (pdb_content and selection):
        return None

    chain_selection: dict[str, str] = {}
    pdb_data: list[dict[str, str]] = extract_chains_from_pdb(file_content=pdb_content)

    for chain_data in pdb_data:
        chain_id: str = chain_data["chain_id"]
        if not chain_id in selection.keys():
            continue

        sequence: str = chain_data["sequence"][selection[chain_id][0]:selection[chain_id][1]]
        chain_selection[chain_id] = sequence

    return chain_selection


def generate_linked_chains(pdb_content: str, chain_selection: dict[str, str], linkage: dict[str, list[str]], linker: str=("GGGGS" * 5)) -> dict[str, str] | None:
    """
    Generates sequences by linking selected residues in a specified order, using a provided linker sequence.
    :param pdb_content: The full pdb file as a string (although not directly used in this function, it's passed from calling functions).
    :param chain_selection: A dictionary where keys are chain IDs and values are the selected sequence strings for those chains.
    :param linkage: A dictionary defining the order of chains to be linked. Keys represent the final linked chain IDs, and values are lists of original chain IDs to be concatenated.
    :param linker: The amino acid sequence to be used as a linker between different chain selections (default is "GGGGS" repeated 5 times).
    :return: A dictionary of linked chain amino acid sequences, or None if input is invalid or a specified chain in linkage is not found.
    """
    if not (pdb_content and chain_selection and linkage):
        return None

    chains: dict[str, str] = {key: "" for key in linkage.keys()}

    for mesa_chain_id in chains:
        for i, chain_id in enumerate(linkage[mesa_chain_id]):
            if chain_id in chain_selection:
                chains[mesa_chain_id] += chain_selection[chain_id]

                if i < len(linkage[mesa_chain_id]) - 1:
                    chains[mesa_chain_id] += linker

            else:
                print(f"Warning: Chain '{chain_id}' specified in linkage for '{mesa_chain_id}' not found in chain_selection.")
                return None

    return chains