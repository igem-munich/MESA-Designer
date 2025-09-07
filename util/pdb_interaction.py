from Bio import SeqIO
from io import StringIO
import re
import requests

def extract_chains_from_pdb(file_path: str | None = None, file_content: str | None = None) -> list[dict]:
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
