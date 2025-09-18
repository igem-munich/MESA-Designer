from __future__ import annotations
import copy
from typing import Dict, List, Optional, Any
from .part import AnnotatedPart, Annotation
from mesa_designer import TMD_DATA, AIP_DATA, SIGNAL_SEQS, NTEV_DATA, CTEV_DATA, TEVP_DATA, PRS_DATA, FRET_ICDs

# Defines the canonical order of components within a MESA chain for assembly and annotation purposes.
MESA_ORDER: List[str] = ["signal_peptide",
                         "tags",
                         "binder",
                         "tmd_linker",
                         "tmd",
                         "fret",
                         "protease",
                         "prs",
                         "cargo",
                         "aip"]


class MesaChain:
    """
    Represents a single MESA (Modularized Extracellular Sensing Assembly) chain,
    composed of various biological parts in a defined order.
    """
    def __init__(self, name: Optional[str] = None) -> None:
        """
        Initializes a new MesaChain instance.
        :param name: An optional name for the MESA chain.
        :return: None
        """
        self.name: Optional[str] = name
        # A dictionary to store the individual AnnotatedPart components of the chain,
        # where keys are component names (e.g., "binder", "tmd") and values are AnnotatedPart objects.
        self.parts: Dict[str, AnnotatedPart] = {}

    def get_parts(self) -> Dict[str, AnnotatedPart]:
        """
        Retrieves the dictionary of parts currently comprising this MESA chain.
        :return: A dictionary where keys are component names and values are AnnotatedPart objects.
        """
        return self.parts

    def add_binder(self, sequence: str, name: Optional[str] = None, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a binder component to the MESA chain.
        :param sequence: The amino acid sequence of the binder.
        :param name: An optional name for the binder part. Defaults to "Binder".
        :param annotation: An optional annotation string for the binder. Defaults to "Binder".
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the sequence is empty or None.
        """
        if not sequence:
            raise ValueError("Sequence cannot be None")

        sequence = sequence.upper()

        self.parts["binder"] = AnnotatedPart(sequence=sequence,
                                             name=name if name else "Binder",
                                             seq_annotations=[Annotation("Binder" if not annotation else annotation, sequence=sequence)])

        return self

    def add_tmd_linker(self, sequence: Optional[str] = None, name: Optional[str] = None, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a TMD (Transmembrane Domain) linker component to the MESA chain.
        :param sequence: The amino acid sequence of the linker. If None, defaults to "GGGS" * 10.
        :param name: An optional name for the linker part. Defaults to "TMD Linker".
        :param annotation: An optional annotation string for the linker. Defaults to "Linker".
        :return: The MesaChain instance, allowing for method chaining.
        """
        if not sequence:
            sequence = "GGGS" * 10

        sequence = sequence.upper()

        self.parts["tmd_linker"] = AnnotatedPart(sequence=sequence,
                                                 name=name if name else "TMD Linker",
                                                 seq_annotations=[Annotation("Linker" if not annotation else annotation, sequence=sequence)])

        return self

    def add_protease(self, protease_name: str) -> MesaChain:
        """
        Adds a predefined protease component to the MESA chain.
        :param protease_name: The name of the protease to add (must be in NTEV_DATA, CTEV_DATA, or TEVP_DATA).
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the protease_name is not found in the predefined protease data.
        """
        proteases: Dict[str, Any] = NTEV_DATA
        proteases.update(CTEV_DATA)
        proteases.update(TEVP_DATA)

        if protease_name not in proteases.keys():
            raise ValueError(
                f"{protease_name} is not a valid Protease name. Please only use available NTEV, CTEV or TEVP names")

        # The sequence is assumed to be the second element (index 1) in the protease data tuple/list.
        self.parts["protease"] = AnnotatedPart(sequence=proteases[protease_name][1],
                                               name=protease_name,
                                               seq_annotations=[Annotation(name=protease_name, sequence=proteases[protease_name][1])])

        return self

    def add_custom_protease(self, sequence: str, name: Optional[str] = None, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a custom protease component to the MESA chain using a provided sequence.
        :param sequence: The amino acid sequence of the custom protease.
        :param name: An optional name for the custom protease part. Defaults to "Protease".
        :param annotation: An optional annotation string for the custom protease. Defaults to "Protease".
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the sequence is empty or None.
        """
        if not sequence:
            raise ValueError("Sequence cannot be None")

        sequence = sequence.upper()

        self.parts["protease"] = AnnotatedPart(sequence=sequence,
                                               name=name if name else "Protease",
                                               seq_annotations=[Annotation("Protease" if not annotation else annotation, sequence=sequence)])

        return self

    def add_prs(self, prs_name: Optional[str] = None) -> MesaChain:
        """
        Adds a predefined Protease Recognition Site (PRS) component to the MESA chain.
        :param prs_name: The name of the PRS to add. If None, defaults to "PRS". Must be found in PRS_DATA.
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the prs_name is not found in the predefined PRS data.
        """
        if not prs_name:
            prs_name = "PRS"

        if prs_name not in PRS_DATA.keys():
            raise ValueError(f"{prs_name} is not a valid PRS name. Please only use available PRS names")

        # The sequence is assumed to be the second element (index 1) in the PRS data tuple/list.
        self.parts["prs"] = AnnotatedPart(sequence=PRS_DATA[prs_name][1],
                                          name=prs_name,
                                          seq_annotations=[Annotation(prs_name, sequence=PRS_DATA[prs_name][1])])

        return self

    def add_custom_prs(self, sequence: str, name: Optional[str] = None, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a custom Protease Recognition Site (PRS) component to the MESA chain.
        :param sequence: The amino acid sequence of the custom PRS.
        :param name: An optional name for the custom PRS part. Defaults to "PRS".
        :param annotation: An optional annotation string for the custom PRS. Defaults to "PRS".
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the sequence is empty or None.
        """
        if not sequence:
            raise ValueError("Protease Recognition Sequence cannot be None")

        sequence = sequence.upper()

        self.parts["prs"] = AnnotatedPart(sequence=sequence,
                                          name=name if name else "PRS",
                                          seq_annotations=[Annotation(name="PRS" if not annotation else annotation, sequence=sequence)])

        return self

    def add_cargo(self, sequence: str, name: Optional[str] = None, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a cargo (target) component to the MESA chain.
        :param sequence: The amino acid sequence of the cargo.
        :param name: An optional name for the cargo part. Defaults to "Cargo".
        :param annotation: An optional annotation string for the cargo. Defaults to "Cargo".
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the sequence is empty or None.
        """
        if not sequence:
            raise ValueError("Cargo Sequence cannot be None")

        sequence = sequence.upper()

        self.parts["cargo"] = AnnotatedPart(sequence=sequence,
                                            name=name if name else "Cargo",
                                            seq_annotations=[Annotation(name="Cargo" if not annotation else annotation, sequence=sequence)])

        return self

    def add_tmd(self, tmd_name: str) -> MesaChain:
        """
        Adds a predefined Transmembrane Domain (TMD) component to the MESA chain.
        :param tmd_name: The name of the TMD to add. Case-insensitive, must be found in TMD_DATA.
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the tmd_name is not found in the predefined TMD data.
        """
        tmd_name = tmd_name.upper()
        if tmd_name not in TMD_DATA.keys():
            raise ValueError(
                f"{tmd_name} is not a valid TMD name. Please only use available TMD names or use a custom TMD")

        # The sequence is assumed to be the second element (index 1) in the TMD data tuple/list.
        self.parts["tmd"] = AnnotatedPart(sequence=TMD_DATA[tmd_name][1],
                                          name=f"{tmd_name}_TMD",
                                          seq_annotations=[Annotation(f"{tmd_name}_TMD", sequence=TMD_DATA[tmd_name][1])])

        return self

    def add_custom_tmd(self, sequence: str, name: Optional[str] = None, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a custom Transmembrane Domain (TMD) component to the MESA chain.
        :param sequence: The amino acid sequence of the custom TMD.
        :param name: An optional name for the custom TMD part. Defaults to "TMD".
        :param annotation: An optional annotation string for the custom TMD. Defaults to "{name}_TMD".
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the sequence is empty or None.
        """
        if not sequence:
            raise ValueError("TMD Sequence cannot be None")

        sequence = sequence.upper()

        self.parts["tmd"] = AnnotatedPart(sequence=sequence,
                                          name=name if name else "TMD",
                                          seq_annotations=[Annotation(f"{name}_TMD" if not annotation else annotation, sequence=sequence)])

        return self

    def add_signal_peptide(self, peptide_name: Optional[str] = None) -> MesaChain:
        """
        Adds a predefined signal peptide component to the MESA chain.
        :param peptide_name: The name of the signal peptide to add. If None, defaults to "CD4". Case-insensitive, must be found in SIGNAL_SEQS.
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the peptide_name is not found in the predefined signal sequence data.
        """
        if not peptide_name:
            peptide_name = "CD4"

        peptide_name = peptide_name.upper()
        if peptide_name not in SIGNAL_SEQS.keys():
            raise ValueError(
                f"{peptide_name} is not a valid peptide name. Please only use available peptide names or use a custom peptide sequence")

        # The sequence is assumed to be the second element (index 1) in the signal sequence data tuple/list.
        self.parts["signal_peptide"] = AnnotatedPart(sequence=SIGNAL_SEQS[peptide_name][1],
                                                     name=f"{peptide_name}_Signal_Peptide",
                                                     seq_annotations=[Annotation(f"{peptide_name}_Signal_Peptide", sequence=SIGNAL_SEQS[peptide_name][1])])

        return self

    def add_custom_signal_peptide(self, sequence: str, name: Optional[str] = None, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a custom signal peptide component to the MESA chain.
        :param sequence: The amino acid sequence of the custom signal peptide.
        :param name: An optional name for the custom signal peptide part. Defaults to "Signal_Peptide".
        :param annotation: An optional annotation string for the custom signal peptide. Defaults to "{annotation}_Signal_Peptide".
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the sequence is empty or None.
        """
        if not sequence:
            raise ValueError("Peptide sequence cannot be None")

        sequence = sequence.upper()

        self.parts["signal_peptide"] = AnnotatedPart(sequence=sequence,
                                                     name=name if name else "Signal_Peptide",
                                                     seq_annotations=[Annotation(f"{annotation}_Signal_Peptide" if not annotation else annotation, sequence=sequence)])

        return self

    def add_aip(self, aip_name: str) -> MesaChain:
        """
        Adds a predefined Autoinducing Peptide (AIP) component to the MESA chain.
        :param aip_name: The name of the AIP to add. Case-insensitive, must be found in AIP_DATA.
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the aip_name is not found in the predefined AIP data.
        """
        aip_name = aip_name.upper()
        if aip_name not in AIP_DATA.keys():
            raise ValueError(
                f"{aip_name} is not a valid AIP name. Please only use available AIP names or use a custom AIP")

        # The sequence is assumed to be the second element (index 1) in the AIP data tuple/list.
        self.parts["aip"] = AnnotatedPart(sequence=AIP_DATA[aip_name][1],
                                          name=f"{aip_name}_AIP",
                                          seq_annotations=[Annotation(f"{aip_name}_AIP", sequence=AIP_DATA[aip_name][1])])

        return self

    def add_custom_aip(self, sequence: str, name: str, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a custom Autoinducing Peptide (AIP) component to the MESA chain.
        :param sequence: The amino acid sequence of the custom AIP.
        :param name: The name for the custom AIP part.
        :param annotation: An optional annotation string for the custom AIP. Defaults to "{name}_AIP".
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the sequence is empty or None.
        """
        if not sequence:
            raise ValueError("AIP Sequence cannot be None")

        sequence = sequence.upper()

        self.parts["aip"] = AnnotatedPart(sequence=sequence,
                                          name=name,
                                          seq_annotations=[Annotation(f"{name}_AIP" if not annotation else annotation, sequence=sequence)])

        return self

    def add_component(self, name: str, sequence: str, part_id: Optional[str] = None, description: Optional[str] = None, annotation: Optional[str] = None) -> MesaChain:
        """
        Adds a generic custom component to the MESA chain.
        :param name: The name of the component, which will also be used as the key in the parts dictionary.
        :param sequence: The amino acid sequence of the component.
        :param part_id: An optional ID for the component part. Defaults to 'name'.
        :param description: An optional description for the component part. Defaults to 'name'.
        :param annotation: An optional annotation string for the component. Defaults to 'name'.
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the name or sequence is empty or None.
        """
        if not name:
            raise ValueError("Name cannot be None")

        if not sequence:
            raise ValueError("Sequence cannot be None")

        sequence = sequence.upper()

        self.parts[name] = AnnotatedPart(sequence=sequence,
                                                          name=name,
                                                          part_id=part_id if part_id else name,
                                                          description=description if description else name,
                                                          seq_annotations=[Annotation(name if not annotation else annotation, sequence=sequence)])

        return self

    def add_part(self, name: str, annotated_part: AnnotatedPart) -> MesaChain:
        """
        Adds an already constructed AnnotatedPart object to the MESA chain under a specified name.
        :param name: The name to assign to the AnnotatedPart within the MESA chain's parts dictionary.
        :param annotated_part: The AnnotatedPart object to add.
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the name or annotated_part is None.
        """
        if name is None:
            raise ValueError("Name cannot be None")

        if annotated_part is None:
            raise ValueError("Annotated Part cannot be None")

        self.parts[name] = annotated_part

        return self

    def remove_component(self, component: str) -> MesaChain:
        """
        Removes a component from the MESA chain by its name.
        :param component: The name of the component to remove.
        :return: The MesaChain instance, allowing for method chaining.
        :raises ValueError: If the component name is empty or None, or if the component is not in the chain.
        """
        if not component:
            raise ValueError("Component cannot be None")

        if component in self.parts.keys():
            self.parts.pop(component)

        else:
            raise ValueError("Component not in current MESA Chain")

        return self

    def to_annotated_part(self, name: str, part_id: Optional[str] = None, description: Optional[str] = None) -> AnnotatedPart:
        """
        Converts the MesaChain into a single AnnotatedPart, concatenating all its components
        in the order defined by MESA_ORDER and adding appropriate linkers.
        :param name: The name for the resulting AnnotatedPart.
        :param part_id: An optional ID for the resulting AnnotatedPart. Defaults to 'name'.
        :param description: An optional description for the resulting AnnotatedPart. Defaults to 'name'.
        :return: A single AnnotatedPart representing the entire MESA chain.
        """
        temp: MesaChain = copy.deepcopy(self)

        sequence: str = ""
        seq_annotations: List[Annotation] = []
        for component in MESA_ORDER:
            if component in temp.parts.keys():
                # Shift annotations of the current part to their correct position in the concatenated sequence.
                # get_annotations()[0] assumes each AnnotatedPart in self.parts has at least one annotation.
                seq_annotations.append(temp.parts[component].get_annotations()[0].shift_annotation(len(sequence)))
                # Concatenate the sequence of the current part.
                sequence += temp.parts[component].get_sequence()

                # Conditionally add short linker sequences (GGGSGGGS) between specific components.
                match component:
                    case "binder":
                        # Add linker after binder if no explicit TMD linker is present.
                        if "tmd_linker" not in temp.parts.keys():
                            sequence += "GGGSGGGS"

                    case "tmd":
                        # Add linker after TMD.
                        sequence += "GGGSGGGS"

                    case "protease":
                        # Add linker after protease if PRS, cargo, or AIP are present.
                        if "prs" in temp.parts.keys() or "cargo" in temp.parts.keys() or "aip" in temp.parts.keys():
                            sequence += "GGGSGGGS"

                    case "prs":
                        # Add linker after PRS if cargo is present.
                        if "cargo" in temp.parts.keys():
                            sequence += "GGGSGGGS"

                    case "cargo":
                        # Add linker after cargo if AIP is present.
                        if "aip" in temp.parts.keys():
                            sequence += "GGGSGGGS"

                        # Add a stop codon (*) after cargo.
                        sequence += "*"

        # Ensure the final sequence starts with a Methionine (M) and adjust annotations if necessary.
        if not sequence.startswith("M"):
            sequence = "M" + sequence
            for annotation in seq_annotations:
                annotation.shift_annotation(1)  # Shift all existing annotations by 1 due to added 'M'.

        sequence = sequence.upper()

        return AnnotatedPart(sequence=sequence,
                             name=name,
                             part_id=part_id if part_id else name,
                             description=description if description else name,
                             seq_annotations=seq_annotations)

    def to_genbank_string(self) -> str:
        """
        Converts the MesaChain into a GenBank formatted string.
        :return: A string representation of the MesaChain in GenBank format.
        """
        # Converts the chain to an AnnotatedPart and then calls its to_genbank_string method.
        return self.to_annotated_part(name="mesa_chain").to_genbank_string()

    def save_genbank_file(self, file_path: str, name: Optional[str] = None) -> None:
        """
        Saves the MesaChain as a GenBank file.
        :param file_path: The full path where the GenBank file should be saved.
        :param name: An optional name for the construct within the GenBank file. If None, the name is derived from the file_path.
        :return: None
        """
        # Derives the name for the AnnotatedPart from the file path if not provided.
        # Then converts to an AnnotatedPart and saves it.
        self.to_annotated_part(name=name if name else file_path.split("/")[-1].split(".")[-2]).save_genbank_file(file_path)

    def to_fret_chains(self) -> 'MesaAssembly':
        """
        Transforms the current MesaChain into a MesaAssembly containing two FRET-compatible
        chains: one with mVenus and one with mCerulean, replacing certain components.
        Protease, cargo, PRS, and AIP components are removed from the base chain before
        adding FRET components.
        :return: A MesaAssembly containing the two FRET-modified MesaChains.
        """
        chains: MesaAssembly = MesaAssembly()
        chain: MesaChain = copy.deepcopy(self)
        # Remove components that are usually cleaved off or not relevant for FRET constructs.
        for component in {"protease", "cargo", "prs", "aip"}:
            if component in chain.parts.keys():
                chain.remove_component(component)

        # Create the mVenus FRET chain.
        chains.set_chain("mVenus", copy.deepcopy(chain).add_component("fret", sequence=FRET_ICDs["mVenus"][1], annotation="mVenus"))
        # Create the mCerulean FRET chain (using the modified base chain, then adding mCerulean).
        chains.set_chain("mCerulean", chain.add_component("fret", sequence=FRET_ICDs["mCerulean"][1], annotation="mCerulean"))

        return chains


class MesaAssembly:
    """
    Represents an assembly of multiple MesaChain objects, typically used for managing
    related constructs, such as FRET pairs.
    """
    def __init__(self, mesa_chains: Optional[Dict[str, MesaChain]] = None) -> None:
        """
        Initializes a new MesaAssembly instance.
        :param mesa_chains: An optional dictionary of MesaChain objects, where keys are names for the chains and values are MesaChain instances.
        :return: None
        """
        self.mesa_chains: Dict[str, MesaChain] = mesa_chains if mesa_chains else {}

    def set_chain(self, name: str, mesa_chain: MesaChain) -> 'MesaAssembly':
        """
        Adds or updates a MesaChain within the assembly.
        :param name: The name to assign to the MesaChain in the assembly.
        :param mesa_chain: The MesaChain object to add.
        :return: The MesaAssembly instance, allowing for method chaining.
        :raises ValueError: If the name or mesa_chain is None.
        """
        if not (name and mesa_chain):
            raise ValueError("MESA Chain or Name cannot be None")

        self.mesa_chains[name] = mesa_chain

        return self

    def to_genbank_strings(self) -> Dict[str, str]:
        """
        Converts all MesaChain objects in the assembly into a dictionary of GenBank formatted strings.
        :return: A dictionary where keys are the chain names and values are their GenBank string representations.
        """
        return {name: mesa_chain.to_genbank_string() for name, mesa_chain in self.mesa_chains.items()}

    def save_genbank_files(self, file_path: str) -> None:
        """
        Saves each MesaChain in the assembly as a separate GenBank file.
        The file names will be derived from the chain names within the specified directory.
        :param file_path: The base directory path where the GenBank files should be saved. Each file will be named "{file_path}/{chain_name}.gb".
        :return: None
        """
        for name, chain in self.mesa_chains.items():
            chain.save_genbank_file(file_path + "/" + name + ".gb")
