from fastapi import FastAPI, File, HTTPException, UploadFile, status, Path, Query
from pathlib import Path as pathlibPath
import sys
from pydantic import BaseModel, Field

# add main directory to system path for imports
current_dir = pathlibPath(__file__).resolve().parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

from util.antibody_search import search_antibodies_api
from util.pdb_interaction import get_pdb_from_rcsb, get_fasta_from_rcsb, extract_chains_from_pdb, generate_chain_selection, generate_linked_chains
from util import TMD_DATA, CTEV_DATA, NTEV_DATA, TEVP_DATA, PRS_DATA, SIGNAL_SEQS, PRS_DATA

# pydantic models
class ChainSelection(BaseModel):
    selection: dict[str, tuple[int, int]] = Field(
        ...,
        examples=[
            {"A": (0, 200), "C": (45, 87)},
        ],
        description="Dictionary where keys are chain IDs and values are tuples of 0-indexed start and endpoint of chain selection."
    )


class ChainLinkage(BaseModel):
    linkage: dict[str, list[str]] = Field(
        ...,
        examples=[
            {"A": ["H", "L"], "B": ["B"]},
            {"B": ["L", "H"]}
        ],
        description="Dictionary where A, B represent the chains of a MESA receptor. The values lists representing the order in which chainselections from a provided PDB structure should be linked."
    )


class SplitProteaseAttachmentInput(BaseModel):
    sequences: dict[str, str] = Field(
        ...,
        examples=[{"n": "SEQUENCE_TO_N_TERMINUS", "c": "SEQUENCE_TO_C_TERMINUS"}],
        description="A dictionary with 'n' and 'c' as keys, mapping to the amino acid sequences to which the N and C terminal parts of the protease are to be attached, respectively."
    )
    protease_splits: dict[str, str] = Field(
        {
            "n": "ESLFKGPRDYNPISSTICHLTNESDGHTTSLYGIGFGPFIITNKHLFRRNNGTLLVQSLHGVFKVKNTTTLQQSLIDGRDMIIIRMPKDFPPFPQKLKFREPQREERICLVTTNFQT",
            "c": "KSMSSMVSDTSCTFPSSDGIFWKHWIQTKDGQCGSPLVSTRDGFIVGIHSASNFTNTNNYFTSVPKNFMELKTNQEAQQWVSGWRLNADSVLWGGHKVFMVKPEEPFQPVKEATQLMN"
        },
        description="A dictionary with 'n' and 'c' as keys, mapping to the corresponding split protease amino acid sequences. Defaults to TEV protease splits."
    )

class ProteaseAttachmentInput(BaseModel):
    sequence: str = Field(
        ...,
        examples=["MASEQUENCE"],
        description="The amino acid sequence to attach the protease to."
    )
    protease_sequence: str = Field(
        TEVP_DATA["TEVp"][1],
        examples=[TEVP_DATA["TEVp"][1]],
        description="The protease sequence to attach to the provided sequence. Defaults to the TEV protease sequence."
    )

class PrsAttachmentInput(BaseModel):
    sequences: dict[str, str] = Field(
        ...,
        examples=[{"A": "SEQUENCE_FOR_CHAIN_A", "B": "SEQUENCE_FOR_CHAIN_B"}],
        description="A dictionary mapping identifiers (e.g., chain names) to the amino acid sequences to which the protease recognition sequence will be attached."
    )
    prs_sequence: str = Field(
        PRS_DATA["PRS"][1],
        examples=[PRS_DATA["PRS"][1]],
        description="The Protease Recognition Sequence (PRS) to attach to the provided sequences. Defaults to a standard PRS sequence."
    )



# functions
def get_pdb_with_http_error(pdb_id: str) -> str:
    pdb = get_pdb_from_rcsb(pdb_id)
    if not pdb:
        raise HTTPException(status_code=404, detail=f"PDB ID '{pdb_id}' not found or could not be retrieved from RCSB.")

    return pdb


# init fast api
app = FastAPI(
    title="MESA-Designer API",
    description="API for MESA-Designer to programmatically interact with its functionalities.",
    version="0.1.0",
)


@app.get(path="/", summary="Root endpoint for MESA-Designer API")
async def read_root() -> dict[str, str]:
    """
    Root endpoint that returns a welcome message.
    """
    return {"message": "Welcome to the MESA-Designer API!"}


@app.get(path="/search_antigen", summary="Search for antibodies based on an antigen query")
async def search_antigens(antigen: str = Query(..., description="The antigen name to search for")) -> dict[str, dict[str, str | int | None] | str]:
    """
    Searches the antibody database for entries matching the provided antigen query.
    """
    if not antigen:
        raise HTTPException(status_code=400, detail="Antigen query cannot be empty.")

    results = search_antibodies_api(antigen)
    if not results["sabdab_data"]:
        raise HTTPException(status_code=404, detail=f"No antibodies found for antigen: {antigen}")

    return results


@app.get(path="/pdb/{pdb_id}_chains", summary="Retrieve PDB chain data")
async def get_pdb_chains(pdb_id: str = Path(..., description="The PDB ID to retrieve chain data for")) -> dict[str, str | list[dict[str, str]]]:
    """
    Retrieves the PDB file content for a given PDB ID from RCSB and then extracts and returns detailed information about each chain within the PDB.
    """
    pdb_content: str = get_pdb_with_http_error(pdb_id)

    chains_data = extract_chains_from_pdb(file_content=pdb_content)
    if not chains_data:
        raise HTTPException(status_code=404, detail=f"No chain data extracted for PDB ID '{pdb_id}'.")

    return {"pdb_id": pdb_id, "chains": chains_data}


@app.get(path="/pdb/{pdb_id}_structure", summary="Retrieve PDB file from RCSB")
async def get_pdb_structure(pdb_id: str = Path(..., description="The PDB ID to retrieve PDB file from RCSB for. It should be noted that it is likely faster to get this structure directly from RCSB")) -> dict[str, str]:
    """
    Retrieves a PDB file from RCSB database based on pdb id. It is usually faster to retrieve this from RCSB directly.
    """
    pdb_content: str = get_pdb_with_http_error(pdb_id)

    return {"pdb_id": pdb_id, "pdb_content": pdb_content}


@app.get(path="/pdb/{pdb_id}_fasta", summary="Retrieve FASTA file from RCSB")
async def get_pdb_fasta(pdb_id: str = Path(..., description="The PDB ID to retrieve FASTA file from RCSB for. It should be noted that it is likely faster to get this file directly from RCSB")) -> dict[str, str]:
    """
    Retrieves the FASTA file for a corresponding PDB ID from RCSB database based on pdb id. It is usually faster to retrieve this from RCSB directly.
    """
    pdb_fasta: str | None = get_fasta_from_rcsb(pdb_id)
    if not pdb_fasta:
        raise HTTPException(status_code=404, detail=f"PDB ID '{pdb_id}' not found or could not be retrieved from RCSB.")

    return {"pdb_id": pdb_id, "pdb_content": pdb_fasta}


@app.post(path="/pdb/{pdb_id}/generate_chain_selection", summary="Generate chain selection from selected PDB cahins and residues")
async def get_chain_selection(selection_data: ChainSelection,
                              pdb_id: str = Path(..., description="The PDB ID for which to generate the chain selection")) -> dict[str, str | dict[str, str]]:
    """
    Creates a dictionary representing the chain selection in the specified pdb structure. This can be used in later requests to specify selected chains
    """
    pdb_content: str = get_pdb_with_http_error(pdb_id)

    chain_selection: dict[str, str] | None = generate_chain_selection(pdb_content, selection_data.selection)
    if not chain_selection:
        raise HTTPException(status_code=400, detail=f"Could not select chains from pdb with the given selection data. Ensure valid chain IDs and residue numbers.")

    return {"pdb_id": pdb_id, "chain_selection": chain_selection}


@app.post(path="/pdb/{pdb_id}/generate_linked_chains", summary="Link selected chains/residues from a specified pdb structure in supplied order and attach custom linkers")
async def get_linked_chains(selection_data: ChainSelection,
                            linkage_data: ChainLinkage,
                            pdb_id: str = Path(..., description="The PDB ID which supplies the selected chains"),
                            linker: str = Query(("GGGGS" * 5), description="The amino acid sequence to use as a linker between concatenated chains. Defaults to 5 repeats of 'GGGGS'.")
                            ) -> dict[str, str | dict[str, str] | dict[str, list[str]]]:
    """
    Generates sequences resulting from linking selected residues in specified order
    """
    pdb_content: str = get_pdb_with_http_error(pdb_id)

    chain_selection: dict[str, str] | None = generate_chain_selection(pdb_content, selection_data.selection)
    if not chain_selection:
        raise HTTPException(status_code=400, detail=f"Could not select chains from pdb with the given selection data. Ensure valid chain IDs and residue numbers.")

    chains: dict[str, str] | None = generate_linked_chains(pdb_content, chain_selection, linkage_data.linkage, linker)
    if not chains:
        raise HTTPException(status_code=400, detail=f"Could not link chains from pdb with given selection and linkage data. Ensure valid chain IDs and residue numbers.")

    return {"pdb_id": pdb_id,
            "chain_selection": chain_selection,
            "linkage_data": linkage_data.linkage,
            "linker": linker,
            "mesa_chains": chains}

@app.get(path="/tmd/overview", summary="Receive an overview over a selection of TMDs taken from various papers.")
async def get_tmd_overview() -> dict[str, dict[str, list[str]] | list[str]]:
    """
    Receive an overview over a selection of common TMDs taken from various papers and associated paper references.
    """
    return {
        "tmd_data": TMD_DATA,
        "sources": ["https://academic.oup.com/synbio/article/5/1/ysaa017/5913400", "https://pubs.acs.org/doi/10.1021/sb400128g"]}


@app.post(path="/tmd/attach", summary="Attach a specified TMD to a provided chain with specified linker and prepend CD4 signal peptide.")
async def get_attached_tmd(tmd: str = Query(..., description="The ID of the TMD to attach (e.g., 'FGFR4')."),
                           sequence: str = Query(..., description="The amino acid sequences to which the TMD should be attached. Usually the result of a post request to /pdb/{pdb_id}/generate_linked_chains."),
                           linker: str = Query(("GGGS" * 10), description="The linker to be used when attaching the TMD sequence.")) -> dict[str, str]:
    """
    Attach specified TMD to provided sequence using specified or default linker.
    """
    if not tmd.upper() in [key.upper() for key in TMD_DATA.keys()]:
        raise HTTPException(status_code=400, detail="Invalid TMD provided. You an get an overview of available options at /tmd/overview")

    tmd_sequence = TMD_DATA[tmd.upper()][1].upper()
    return {"cd4": SIGNAL_SEQS["CD4"][1], "tmd": tmd.upper(), "tmd_sequence": tmd_sequence, "linker": linker.upper(), "combined_sequence": SIGNAL_SEQS["CD4"][1] + sequence.upper() + linker.upper() + tmd_sequence}


@app.get(path="/protease/overview", summary="Get an overview over protease design options and select amino acid sequences")
async def get_protease_overview() -> dict[str, list[str] | dict[str, str | dict[str, str | list[str] | dict[str, list[str]]]]]:
    """
    Return an overview of available protease designs and select amino acid sequences. Additionally, point to papers which outline these and link to tools for custom protease design.
    """
    return {
        "separate_chains": {
            "description": "This design uses one mesa chain as the protease chain and the other chain as the cargo chain.",
            "TEVp_sequences": TEVP_DATA
        },
        "split_protease": {
            "description": "This design uses a split protease. This means that the protease is split and the separate parts are attached to separate chains. These reconstitute upon binder dimerization.",
            "NTEVp_sequences": NTEV_DATA,
            "CTEVp_sequences": CTEV_DATA,
            "spell_tool": {
                "description": "This tool can be used to generally guide the splitting process of proteins and enzymes.",
                "link": "https://dokhlab.med.psu.edu/spell/login.php"
            }
        },
        "sources": ["https://academic.oup.com/synbio/article/5/1/ysaa017/5913400", "https://pubs.acs.org/doi/10.1021/sb400128g"],
        "common": {
            "prs": {
                "name": "protease recognition sequence",
                "description": "The protease recognition sequence is the amino acid sequence which is cleaved by the protease (split or whole) upon dimerization.",
                "prs_sequences": PRS_DATA
            },
            "aip": {
                "name": "auto-inhibitory peptide",
                "description": "The AIP sequence can reduce background signaling by reversibly blocking the active site of the protease"
            }
        }
    }


@app.get(path="/protease/split_tev_protease")
async def get_split_tevp() -> dict[str, dict[str, list[str]]]:
    return {
        "ntevp": NTEV_DATA,
        "ctevp": CTEV_DATA
    }


@app.post(path="/protease/attach_split", summary="Attach the chains of a split protease to two uploaded sequences with short linker.")
async def get_attached_split_protease(split_protease_data: SplitProteaseAttachmentInput) -> dict[str, str | dict[str, str]]:
    """
    Generate combined chain sequences of input sequences (commonly received from /tmd/attach or /pdb/{pdb_id}/generate_linked_chains
    """
    sequences = split_protease_data.sequences
    protease_splits = split_protease_data.protease_splits

    if not sequences or any(len(sequence) < 1 for sequence in sequences.values()):
        raise HTTPException(status_code=400, detail="You need to provide chains to attach the split protease components to.")

    return {"protease_splits": protease_splits,
            "n_sequence": sequences["n"] + "GGGSGGGS" + protease_splits["n"],
            "c_sequence": sequences["c"] + "GGGSGGGS" + protease_splits["c"]}


@app.post(path="/protease/attach_protease", summary="Attach a protease sequence to a provided sequence with a short linker.")
async def get_attached_protease(protease_attachment_data: ProteaseAttachmentInput) -> dict[str, str]:
    """
    Generate combined sequence of provided sequence and default TEV protease or provided custom protease
    """
    sequence = protease_attachment_data.sequence
    protease_sequence = protease_attachment_data.protease_sequence

    if not sequence or len(sequence) < 1:
        raise HTTPException(status_code=400, detail="You need to provide a sequence to attach the protease sequence to.")

    return {"protease_sequence": protease_sequence, "combined_sequence": sequence + "GGGSGGGS" + protease_sequence}


@app.post(path="/protease/attach_prs", summary="Attach a protease recognition sequence (PRS) to a provided set of sequences with a short linker.")
async def get_attached_prs(prs_attachment_data: PrsAttachmentInput) -> dict[str, str | dict[str, str]]:
    """
    Generate combined sequences of provided sequences and provided PRS sequence
    """
    sequences = prs_attachment_data.sequences
    prs_sequence = prs_attachment_data.prs_sequence

    if not sequences or any(len(sequence) < 1 for sequence in sequences):
        raise HTTPException(status_code=400, detail="Provided sequences need to be at least of length 1.")

    return {"prs": prs_sequence, "sequences": {chain_id: sequence + "GGGSGGGS" + prs_sequence for chain_id, sequence in sequences.items()}}