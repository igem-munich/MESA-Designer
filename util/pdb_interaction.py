from Bio import SeqIO
from io import StringIO
import re
import requests

def extract_chains_from_pdb(file_path: str | None = None, file_content: str | None = None) -> list[dict[str, str]]:
    chains_data = []
    try:
        records = SeqIO.parse(StringIO(file_content) if file_content else file_path, "pdb-atom")
        for record in records:
            if ":" in record.id:
                chain_id_display = record.id.split(":")[-1]
            else:
                if record.id and record.id[-1].isalpha():
                    chain_id_display = record.id[-1]
                else:
                    chain_id_display = record.id

            fasta_name = re.sub("[^a-zA-Z0-9]", "", record.id.replace(":", "_"))

            chains_data.append({
                "id": record.id,
                "chain_id": chain_id_display,
                "fasta_name": fasta_name,
                "sequence": str(record.seq)
            })
    except Exception as e:
        pass

    return chains_data


def extract_chains_from_fasta(fasta: str) -> list[dict] | None:
    chains_data = []
    try:
        fasta_io = StringIO(fasta)
        records = SeqIO.parse(fasta_io, "fasta")

        for record in records:
            chains_data.append(record)

        return chains_data

    except Exception as e:
        print(e)
        return None


def convert_chains_to_fasta_string(chains_data: list[dict]) -> str:
    """Converts a list of chain dictionaries to a FASTA formatted string."""
    fasta_string = StringIO()
    for chain in chains_data:
        fasta_string.write(f">{chain['fasta_name']}\n{chain['sequence']}\n")
    return fasta_string.getvalue()


def extract_fasta_from_pdb(file_path: str) -> str:
    """
    Original function, modified to use extract_chains_from_pdb for consistency.
    """
    chains_data = extract_chains_from_pdb(file_path)
    return convert_chains_to_fasta_string(chains_data)


def cut_protein(file_content: str) -> str:
    raise Exception("Feature currently under construction!")


def get_pdb_from_rcsb(pdb_id: str) -> str | None:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        res = requests.get(url)
        res.raise_for_status()
        return res.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching PDB file for {pdb_id}: {e}")
        return None


def get_fasta_from_rcsb(pdb_id: str) -> str | None:
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    try:
        res = requests.get(url)
        res.raise_for_status()
        return res.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching FASTA file for {pdb_id}: {e}")
        return None


def generate_chain_selection(pdb_content: str, selection: dict[str, tuple[int, int]]) -> dict[str, str] | None:
    """
    Generates a FASTA string from selected chains and residues within PDB content.
    :param pdb_content: The full pdb file as a string
    :param selection: A tuple of two integers representing the starting residue index (0-indexed, included) and the final residue index (excluded)
    :return: A FASTA string of selected residues and chains
    """
    if not (pdb_content and selection):
        return None

    chain_selection: dict[str, str] = {}
    pdb_data = extract_chains_from_pdb(file_content=pdb_content)

    for chain_data in pdb_data:
        chain_id = chain_data["chain_id"]
        if not chain_id in selection.keys():
            continue

        sequence = chain_data["sequence"][selection[chain_id][0]:selection[chain_id][1]]
        chain_selection[chain_id] = sequence

    return chain_selection


def generate_linked_chains(pdb_content: str, chain_selection: dict[str, str], linkage: dict[str, list[str]], linker: str=("GGGGS" * 5)) -> dict[str, str] | None:
    """
    Generates sequences resulting from linking selected residues in specified order.
    :param pdb_content: The full pdb file as a string.
    :param chain_selection: A dictionary where keys are chain IDs and values are the selected sequence strings for those chains.
    :param linkage: A dictionary of chains representing in what order selections should be linked
    :param linker: The amino acid sequence which should be used as a linker between different chain selections
    :return: A dictionary of A and B chain amino acids after linkage
    """
    if not (pdb_content and chain_selection and linkage):
        return None

    chains = {key: "" for key in linkage.keys()}

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
