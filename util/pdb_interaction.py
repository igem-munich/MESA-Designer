from Bio import SeqIO
import re

def extract_fasta_from_pdb(file_path: str) -> str:
    fasta: str = ""

    records = SeqIO.parse(file_path, "pdb-atom")
    for record in records:
        name = re.sub("[^a-zA-Z0-9]", "", record.id)
        seq = str(record.seq)

        fasta += f">{name}\n{seq}"

    return fasta