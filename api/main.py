from fastapi import FastAPI, HTTPException, Path, Query
from pathlib import Path as pathlibPath
import sys
from pydantic import BaseModel, Field

# Add the main directory of the project to the system path.
# This allows for importing utility modules from the 'util' package.
current_dir = pathlibPath(__file__).resolve().parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

# Import utility functions for antibody searching and PDB interaction from the 'util' package.
from util.antibody_search import search_antibodies_api
from util.pdb_interaction import get_pdb_from_rcsb, get_fasta_from_rcsb, extract_chains_from_pdb, generate_chain_selection, generate_linked_chains
# Import data dictionaries for various biological components like TMDs, proteases, signal sequences, etc.
from util import TMD_DATA, CTEV_DATA, NTEV_DATA, TEVP_DATA, SIGNAL_SEQS, PRS_DATA, AIP_DATA, FRET_ICDs, TAG_SEQS

# Pydantic models are used for data validation and serialization of request bodies and responses.

class ChainSelection(BaseModel):
    """
    Pydantic model to define the structure for selecting specific chain segments
    within a PDB structure.
    """
    selection: dict[str, tuple[int, int]] = Field(
        ...,
        examples=[
            {"A": (0, 200), "C": (45, 87)},
        ],
        description="Dictionary where keys are chain IDs and values are tuples of 0-indexed start and endpoint of chain selection."
    )


class ChainLinkage(BaseModel):
    """
    Pydantic model to define how selected chains from a PDB should be linked together
    to form new MESA receptor chains (Chain A, Chain B).
    """
    linkage: dict[str, list[str]] = Field(
        ...,
        examples=[
            {"A": ["H", "L"], "B": ["B"]},
            {"B": ["L", "H"]}
        ],
        description="Dictionary where A, B represent the chains of a MESA receptor. The values lists representing the order in which chainselections from a provided PDB structure should be linked."
    )


class SplitProteaseAttachmentInput(BaseModel):
    """
    Pydantic model for attaching split protease components to sequences.
    """
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
    """
    Pydantic model for attaching a complete protease sequence to a single provided sequence.
    """
    sequence: str = Field(
        ...,
        examples=["MASEQUENCE"],
        description="The amino acid sequence to attach the protease to."
    )
    protease_sequence: str = Field(
        TEVP_DATA["TEVp"][1], # Default to TEV protease sequence if not provided.
        examples=[TEVP_DATA["TEVp"][1]],
        description="The protease sequence to attach to the provided sequence. Defaults to the TEV protease sequence."
    )


class PrsAttachmentInput(BaseModel):
    """
    Pydantic model for attaching a Protease Recognition Sequence (PRS) to multiple provided sequences.
    """
    sequences: dict[str, str] = Field(
        ...,
        examples=[{"A": "SEQUENCE_FOR_CHAIN_A", "B": "SEQUENCE_FOR_CHAIN_B"}],
        description="A dictionary mapping identifiers (e.g., chain names) to the amino acid sequences to which the protease recognition sequence will be attached."
    )
    prs_sequence: str = Field(
        PRS_DATA["PRS"][1], # Default to a standard TEV-Protease PRS sequence.
        examples=[PRS_DATA["PRS"][1]],
        description="The Protease Recognition Sequence (PRS) to attach to the provided sequences. Defaults to a standard TEV-Protease PRS sequence."
    )


class AipAttachmentInput(BaseModel):
    """
    Pydantic model for attaching an Auto-Inhibitory Peptide (AIP) to multiple provided sequences.
    """
    sequences: dict[str, str] = Field(
        ...,
        examples=[{"A": "SEQUENCE_FOR_CHAIN_A", "B": "SEQUENCE_FOR_CHAIN_B"}],
        description="A dictionary mapping identifiers (e.g., chain names) to the amino acid sequences to which the auto-inhibitory peptide will be attached."
    )
    aip_sequence: str = Field(
        AIP_DATA["AIP"][1], # Default to a standard AIP for TEV-Protease.
        examples=[AIP_DATA["AIP"][1]],
        description="The auto-inhibitory peptide to attach tothe provided sequences. Defaults to a standard AIP for TEV-Protease."
    )


class CargoAttachmentInput(BaseModel):
    """
    Pydantic model for attaching a cargo sequence and optionally a PRS to provided chains.
    """
    sequences: dict[str, str] = Field(
        ...,
        examples=[{"A": "SEQUENCE_FOR_CHAIN_A", "B": "SEQUENCE_FOR_CHAIN_B"}],
        description="A dictionary mapping identifiers (e.g., chain names) to the amino acid sequences to which the cargo will be attached."
    )
    cargo_sequence: str = Field(
        ...,
        examples=[{"AAAAAAAAAAAAAAAAAAAAAA"}],
        description="The desired target's amino acid sequence"
    )
    prepend_prs: bool = Field(
        True,
        description="Optional toggle to prepend a PRS sequence"
    )
    prs_sequence: str = Field(
        PRS_DATA["PRS"][1], # Default PRS sequence if prepended.
        examples=[PRS_DATA["PRS"][1]],
        description="The Protease Recognition Sequence (PRS) to attach to the provided sequences. Defaults to a standard TEV-Protease PRS sequence."
    )


class FretSequenceInput(BaseModel):
    """
    Pydantic model for generating FRET imaging sequences from input sequences.
    """
    sequences: dict[str, str] = Field(
        ...,
        examples=[{"A": "SEQUENCE_FOR_CHAIN_A", "B": "SEQUENCE_FOR_CHAIN_B"}],
        description="A dictionary mapping identifiers (e.g., chain names) to the amino acid sequences to which the cargo will be attached."
    )


class TagSequenceInput(BaseModel):
    """
    Pydantic model for prepending a tag sequence to multiple provided sequences.
    """
    sequences: dict[str, str] = Field(
        ...,
        examples=[{"A": "SEQUENCE_FOR_CHAIN_A", "B": "SEQUENCE_FOR_CHAIN_B"}],
        description="A dictionary mapping identifiers (e.g., chain names) to the amino acid sequences to which the tag will be prepended."
    )
    tag_sequence: str = Field(
        ...,
        description="Amino acid sequence of the tag to be prepended to all supplied sequences"
    )


# Helper functions
def get_pdb_with_http_error(pdb_id: str) -> str:
    """
    Retrieves PDB content from RCSB and raises an HTTPException if retrieval fails.
    :param pdb_id: The ID of the PDB to retrieve.
    :return: The content of the PDB file as a string.
    """
    pdb: str | None = get_pdb_from_rcsb(pdb_id)
    if not pdb:
        raise HTTPException(status_code=404, detail=f"PDB ID '{pdb_id}' not found or could not be retrieved from RCSB.")

    return pdb


# Initialize the FastAPI application.
app: FastAPI = FastAPI(
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
async def search_antigens(antigen: str = Query(..., description="The antigen name to search for")) -> dict[str, list[dict[str, str]] | float]:
    """
    Searches the antibody database for entries matching the provided antigen query.
    :param antigen: The name of the antigen to search for.
    :return: A dictionary containing search results (e.g., SAbDab data) and search duration.
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
    :param pdb_id: The PDB ID for which to retrieve chain data.
    :return: A dictionary containing the PDB ID and a list of dictionaries, each representing a chain with its details.
    """
    pdb_content: str = get_pdb_with_http_error(pdb_id)

    chains_data: list[dict[str, str]] = extract_chains_from_pdb(file_content=pdb_content)
    if not chains_data:
        raise HTTPException(status_code=404, detail=f"No chain data extracted for PDB ID '{pdb_id}'.")

    return {"pdb_id": pdb_id, "chains": chains_data}


@app.get(path="/pdb/{pdb_id}_structure", summary="Retrieve PDB file from RCSB")
async def get_pdb_structure(pdb_id: str = Path(..., description="The PDB ID to retrieve PDB file from RCSB for. It should be noted that it is likely faster to get this structure directly from RCSB")) -> dict[str, str]:
    """
    Retrieves a PDB file from RCSB database based on pdb id. It is usually faster to retrieve this from RCSB directly.
    :param pdb_id: The PDB ID to retrieve the structure for.
    :return: A dictionary containing the PDB ID and its content as a string.
    """
    pdb_content: str = get_pdb_with_http_error(pdb_id)

    return {"pdb_id": pdb_id, "pdb_content": pdb_content}


@app.get(path="/pdb/{pdb_id}_fasta", summary="Retrieve FASTA file from RCSB")
async def get_pdb_fasta(pdb_id: str = Path(..., description="The PDB ID to retrieve FASTA file from RCSB for. It should be noted that it is likely faster to get this file directly from RCSB")) -> dict[str, str]:
    """
    Retrieves the FASTA file for a corresponding PDB ID from RCSB database based on pdb id. It is usually faster to retrieve this from RCSB directly.
    :param pdb_id: The PDB ID to retrieve the FASTA sequence for.
    :return: A dictionary containing the PDB ID and its FASTA content as a string.
    """
    pdb_fasta: str | None = get_fasta_from_rcsb(pdb_id)
    if not pdb_fasta:
        raise HTTPException(status_code=404, detail=f"PDB ID '{pdb_id}' not found or could not be retrieved from RCSB.")

    return {"pdb_id": pdb_id, "pdb_content": pdb_fasta}


@app.post(path="/pdb/{pdb_id}/generate_chain_selection", summary="Generate chain selection from selected PDB cahins and residues")
async def get_chain_selection(selection_data: ChainSelection,
                              pdb_id: str = Path(..., description="The PDB ID for which to generate the chain selection")) -> dict[str, str | dict[str, str]]:
    """
    Creates a dictionary representing the chain selection in the specified pdb structure. This can be used in later requests to specify selected chains.
    :param selection_data: A ChainSelection model specifying the chains and residue ranges.
    :param pdb_id: The PDB ID from which to select chains.
    :return: A dictionary containing the PDB ID and the generated chain selections.
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
    Generates sequences resulting from linking selected residues in specified order.
    :param selection_data: A ChainSelection model specifying the chains and residue ranges.
    :param linkage_data: A ChainLinkage model specifying the order of linking.
    :param pdb_id: The PDB ID supplying the chains.
    :param linker: The amino acid sequence to use as a linker.
    :return: A dictionary containing the PDB ID, chain selection, linkage data, linker, and the assembled MESA chains.
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
    :return: A dictionary containing TMD data and their sources.
    """
    return {
        "tmd_data": TMD_DATA,
        "sources": ["https://academic.oup.com/synbio/article/5/1/ysaa017/5913400", "https://pubs.acs.org/doi/10.1021/sb400128g"]}


@app.post(path="/tmd/attach", summary="Attach a specified TMD to a provided chain with specified linker and prepend CD4 signal peptide.")
async def get_attached_tmd(tmd: str = Query(..., description="The ID of the TMD to attach (e.g., 'FGFR4')."),
                           sequence: str = Query(..., description="The amino acid sequences to which the TMD should be attached. Usually the result of a post request to /pdb/{pdb_id}/generate_linked_chains."),
                           linker: str = Query(("GGGS" * 10), description="The linker to be used when attaching the TMD sequence.")) -> dict[str, str]:
    """
    Attach specified TMD to provided sequence using specified or default linker, and prepend a CD4 signal peptide.
    :param tmd: The ID of the Transmembrane Domain (TMD) to attach.
    :param sequence: The amino acid sequence to which the TMD will be attached.
    :param linker: The amino acid sequence of the linker.
    :return: A dictionary containing the CD4 signal sequence, TMD details, linker, and the combined sequence.
    """
    if not tmd.upper() in [key.upper() for key in TMD_DATA.keys()]:
        raise HTTPException(status_code=400, detail="Invalid TMD provided. You an get an overview of available options at /tmd/overview")

    tmd_sequence: str = TMD_DATA[tmd.upper()][1].upper()
    return {"cd4": SIGNAL_SEQS["CD4"][1], "tmd": tmd.upper(), "tmd_sequence": tmd_sequence, "linker": linker.upper(), "combined_sequence": SIGNAL_SEQS["CD4"][1] + sequence.upper() + linker.upper() + tmd_sequence}


@app.get(path="/protease/overview", summary="Get an overview over protease design options and select amino acid sequences")
async def get_protease_overview() -> dict:
    """
    Return an overview of available protease designs and select amino acid sequences. Additionally, point to papers which outline these and link to tools for custom protease design.
    :return: A dictionary detailing various protease design options, including split and complete proteases,
        and information on associated PRS and AIP sequences, along with relevant sources and tools.
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
                "description": "The AIP sequence can reduce background signaling by reversibly blocking the active site of the protease",
                "aip_sequences": AIP_DATA
            }
        }
    }


@app.get(path="/protease/split_tev_protease")
async def get_split_tevp() -> dict[str, dict[str, list[str]]]:
    """
    Retrieves the sequences for split TEV protease components (N-terminal and C-terminal).
    :return: A dictionary containing the N-terminal (ntevp) and C-terminal (ctevp) TEV protease data.
    """
    return {
        "ntevp": NTEV_DATA,
        "ctevp": CTEV_DATA
    }


@app.post(path="/protease/attach_split_protease", summary="Attach the chains of a split protease to two uploaded sequences with short linker.")
async def get_attached_split_protease(split_protease_data: SplitProteaseAttachmentInput) -> dict[str, str | dict[str, str]]:
    """
    Generate combined chain sequences by attaching the N-terminal and C-terminal parts of a split protease
    to two provided sequences, separated by a short linker ("GGGSGGGS").
    :param split_protease_data: A SplitProteaseAttachmentInput model containing the sequences to attach to and the split protease sequences.
    :return: A dictionary containing the protease splits and the combined N-terminal and C-terminal sequences.
    """
    sequences: dict[str, str] = split_protease_data.sequences
    protease_splits: dict[str, str] = split_protease_data.protease_splits

    if not sequences or any(len(sequence) < 1 for sequence in sequences.values()):
        raise HTTPException(status_code=400, detail="You need to provide chains to attach the split protease components to.")

    # Combine the input sequences with the respective split protease parts and a linker.
    return {"protease_splits": protease_splits,
            "n_sequence": sequences["n"] + "GGGSGGGS" + protease_splits["n"],
            "c_sequence": sequences["c"] + "GGGSGGGS" + protease_splits["c"]}


@app.post(path="/protease/attach_protease", summary="Attach a protease sequence to a provided sequence with a short linker.")
async def get_attached_protease(protease_attachment_data: ProteaseAttachmentInput) -> dict[str, str]:
    """
    Generate a combined sequence by attaching a full protease sequence to a provided sequence,
    separated by a short linker ("GGGSGGGS").
    :param protease_attachment_data: A ProteaseAttachmentInput model containing the sequence to attach to and the protease sequence.
    :return: A dictionary containing the protease sequence and the resulting combined sequence.
    """
    sequence: str = protease_attachment_data.sequence
    protease_sequence: str = protease_attachment_data.protease_sequence

    if not sequence or len(sequence) < 1:
        raise HTTPException(status_code=400, detail="You need to provide a sequence to attach the protease sequence to.")

    # Combine the input sequence with the protease sequence and a linker.
    return {"protease_sequence": protease_sequence, "combined_sequence": sequence + "GGGSGGGS" + protease_sequence}


@app.get(path="/protease/prs_overview", summary="Get an overview over available TEV-Protease recognition sequences.")
async def get_prs_overview() -> dict[str, list[str] | dict[str, list[str]]]:
    """
    Get an overview of pre-selected available TEV-Protease recognition sequences and relevant publications.
    :return: A dictionary containing the PRS sequences and their sources.
    """
    return {"prs_sequences": PRS_DATA, "sources": ["https://academic.oup.com/synbio/article/5/1/ysaa017/5913400", "https://pubs.acs.org/doi/10.1021/sb400128g"]}


@app.post(path="/protease/attach_prs", summary="Attach a protease recognition sequence (PRS) to a provided set of sequences with a short linker.")
async def get_attached_prs(prs_attachment_data: PrsAttachmentInput) -> dict[str, str | dict[str, str]]:
    """
    Generate combined sequences by appending a Protease Recognition Sequence (PRS) to each of the provided sequences,
    separated by a short linker ("GGGSGGGS").
    :param prs_attachment_data: A PrsAttachmentInput model containing the sequences and the PRS sequence to attach.
    :return: A dictionary containing the PRS sequence and a dictionary of resulting combined sequences per chain.
    """
    sequences: dict[str, str] = prs_attachment_data.sequences
    prs_sequence: str = prs_attachment_data.prs_sequence

    if not sequences or any(len(sequence) < 1 for sequence in sequences):
        raise HTTPException(status_code=400, detail="Provided sequences need to be at least of length 1.")

    # For each chain, append the linker and PRS sequence.
    return {"prs_sequence": prs_sequence, "sequences": {chain_id: sequence + "GGGSGGGS" + prs_sequence for chain_id, sequence in sequences.items()}}


@app.get(path="/protease/aip_overview", summary="Get an overview over available TEV-Protease auto-inhibitory peptides and sequences.")
async def get_aip_overview() -> dict[str, list[str] | dict[str, list[str]]]:
    """
    Get an overview of pre-selected available TEV-Protease auto-inhibitory peptides and sequences and relevant publications.
    :return: A dictionary containing AIP sequences and their sources.
    """
    return {"aip_sequences": AIP_DATA, "sources": ["https://academic.oup.com/synbio/article/5/1/ysaa017/5913400", "https://pubs.acs.org/doi/10.1021/sb400128g"]}


@app.post(path="/protease/attach_aip", summary="Attach an auto-inhibitory peptide to a provided set of sequences with a short linker.")
async def get_attached_aip(aip_attachment_data: AipAttachmentInput) -> dict[str, str | dict[str, str]]:
    """
    Generate combined sequences by appending an Auto-Inhibitory Peptide (AIP) to each of the provided sequences,
    separated by a short linker ("GGGSGGGS").
    :param aip_attachment_data: An AipAttachmentInput model containing the sequences and the AIP sequence to attach.
    :return: A dictionary containing the AIP sequence and a dictionary of resulting combined sequences per chain.
    """
    sequences: dict[str, str] = aip_attachment_data.sequences
    aip_sequence: str = aip_attachment_data.aip_sequence

    if not sequences or any(len(sequence) < 1 for sequence in sequences):
        raise HTTPException(status_code=400, detail="Provided sequences need to be at least of length 1.")

    # For each chain, append the linker and AIP sequence.
    return {"aip_sequence": aip_sequence, "sequences": {chain_id: sequence + "GGGSGGGS" + aip_sequence for chain_id, sequence in sequences.items()}}


@app.post(path="/cargo/attach_cargo", summary="Attach a cargo and optionally a protease recognition sequence to provided chains.")
async def get_attached_cargo(cargo_data: CargoAttachmentInput) -> dict[str, bool | str | dict[str, str]]:
    """
    Generate combined sequences by appending a cargo sequence and optionally a Protease Recognition Sequence (PRS)
    to each of the provided sequences, with "GGGSGGGS" linkers in between.
    :param cargo_data: A CargoAttachmentInput model containing the sequences, cargo sequence, and optional PRS details.
    :return: A dictionary containing cargo details, whether PRS was prepended, and the resulting
        combined sequences per chain.
    """
    sequences: dict[str, str] = cargo_data.sequences
    cargo_sequence: str = cargo_data.cargo_sequence
    prepend_prs: bool = cargo_data.prepend_prs
    prs_sequence: str = cargo_data.prs_sequence

    if not sequences or any(len(sequence) < 1 for sequence in sequences):
        raise HTTPException(status_code=400, detail="Provided sequences need to be at least of length 1.")

    if not cargo_sequence or len(cargo_sequence) < 1:
        raise HTTPException(status_code=400, detail="Provided cargo sequence needs to be at least of length 1.")

    if prepend_prs and (not prs_sequence or len(prs_sequence) < 1):
        raise HTTPException(status_code=400, detail="Provided prs sequence needs to be at least of length 1.")

    result: dict[str, str | bool | dict[str, str]] = {"cargo_sequence": cargo_sequence,
                                                      "prepend_prs": prepend_prs,
                                                      "sequences": {
                                                          chain_id: sequence + "GGGSGGGS"
                                                                    + (prs_sequence + "GGGSGGGS" if prepend_prs else "")
                                                                    + cargo_sequence for chain_id, sequence in sequences.items()
                                                      }}

    if prepend_prs:
        result["prs_sequence"] = prs_sequence

    return result


@app.post(path="/extra/fret_sequences", summary="Create FRET imaging sequences from input sequences")
async def get_fret_sequences(fret_data: FretSequenceInput) -> dict[str, dict[str, str]]:
    """
    Generate mVENUS and mCERULEAN FRET imaging sequences for all provided input sequences.
    Each FRET sequence is appended to the input sequence with a "GGGSGGGS" linker.
    :param fret_data: A FretSequenceInput model containing the input sequences.
    :return: A dictionary where the key is "sequences" and the value is a dictionary mapping new chain IDs (e.g., "ChainA_mVenus") to their combined FRET sequences.
    """
    result: dict[str, str] = {}

    # Iterate through each FRET ICD (mVenus, mCerulean) and each provided chain.
    for key in FRET_ICDs.keys():
        for chain_id, sequence in fret_data.sequences.items():
            # Combine the chain sequence with a linker and the FRET ICD sequence.
            result[f"{chain_id}_{key}"] = sequence + "GGGSGGGS" + FRET_ICDs[key][1]

    return {"sequences": result}


@app.get(path="/extra/tag_overview", summary="Get an overview of select Tags from synthetic biology.")
async def get_tag_overview() -> dict[str, list[str] | dict[str, list[str]]]:
    """
    Get an overview of available, commonly used Tag sequences and their sources.
    :return: A dictionary containing tag sequences and their sources.
    """

    return {"sources": ["https://academic.oup.com/synbio/article/5/1/ysaa017/5913400", "https://pubs.acs.org/doi/10.1021/sb400128g"], "tags": TAG_SEQS}


@app.post(path="/extra/attach_tag", summary="Prepend supplied amino acid sequences to all provided sequences.")
async def get_prepended_tag(tag_data: TagSequenceInput) -> dict[str, str | dict[str, str]]:
    """
    Generate new sequences by prepending a specified tag sequence to each of the provided input sequences.
    :param tag_data: A TagSequenceInput model containing the sequences to modify and the tag sequence to prepend.
    :return: A dictionary containing the prepended tag sequence and a dictionary
        of the resulting sequences per chain.
    """
    sequences: dict[str, str] = tag_data.sequences
    tag_sequence: str = tag_data.tag_sequence
    # For each chain, prepend the tag sequence.
    return {"tag_sequence": tag_sequence, "sequences": {chain_id: tag_sequence + sequence for chain_id, sequence in sequences.items()}}