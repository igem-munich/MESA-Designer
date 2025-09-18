from __future__ import annotations
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
from typing import Optional, List


class Annotation:
    """
    Represents an annotation for a biological sequence part,
    detailing its name, type, and genomic coordinates.
    """
    def __init__(self, name: str, sequence: str, start: Optional[int] = None, stop: Optional[int] = None, part_type: Optional[str] = None) -> None:
        """
        Initializes an Annotation object.
        :param name: The name or identifier of the annotation (e.g., "geneX", "CDS_1").
        :param sequence: The biological sequence (as a string) that this annotation applies to. This is used to determine the default 'stop' position if not provided.
        :param start: The start position of the annotation on the sequence (0-indexed, inclusive). If None, defaults to 0.
        :param stop: The stop position of the annotation on the sequence (0-indexed, exclusive). If None, defaults to the length of the provided sequence.
        :param part_type: The type of the annotated feature (e.g., "CDS", "gene", "misc_feature"). If None, defaults to "CDS".
        :return: None
        """
        self.name: str = name
        self.start: int = start if start is not None else 0
        self.stop: int = stop if stop is not None else len(sequence)
        self.type: str = part_type if part_type is not None else "CDS"

    def shift_annotation(self, amount: int) -> Annotation:
        """
        Shifts the start and stop positions of the annotation by a specified amount.
        This is useful for adjusting annotation coordinates when the underlying sequence
        is modified or re-positioned.
        :param amount: The integer value by which to shift the annotation. A positive value shifts it downstream, a negative value upstream.
        :return: The Annotation object itself, with its 'start' and 'stop' attributes updated.
        """
        self.start += amount
        self.stop += amount

        return self

    def to_seq_feature(self) -> SeqFeature:
        """
        Converts the annotation into a Biopython SeqFeature object.
        This allows for easy integration with Biopython's sequence manipulation
        and output functionalities, such as generating GenBank files.
        :return: A Biopython `Bio.SeqFeature.SeqFeature` object representing this annotation. It includes the location, type, and the annotation's name as a qualifier.
        """
        return SeqFeature(location=FeatureLocation(self.start, self.stop), type=self.type, qualifiers={"name": self.name})


class AnnotatedPart:
    """
    Represents a biological sequence part that can have multiple annotations associated with it.
    This class provides methods to manage these annotations and convert the part into
    Biopython SeqRecord objects for GenBank output.
    """
    def __init__(self, sequence: str, name: str, part_id: Optional[str] = None, description: Optional[str] = None, seq_annotations: Optional[List[Annotation]] = None) -> None:
        """
        Initializes an AnnotatedPart object.
        :param sequence: The biological sequence (as a string) of this part.
        :param name: The common name of this part.
        :param part_id: An optional unique identifier for this part. If None, defaults to 'name'.
        :param description: An optional descriptive text for this part. If None, defaults to 'name'.
        :param seq_annotations: An optional list of `Annotation` objects associated with this part. If None, initializes as an empty list.
        :return: None
        """
        self.sequence: str = sequence if sequence else ""
        self.annotations: List[Annotation] = seq_annotations if seq_annotations is not None else []
        self.name: str = name if name else "no_name"
        self.part_id: str = part_id if part_id else name
        self.description: str = description if description else name

    def get_sequence(self) -> str:
        """
        Retrieves the biological sequence string of this part.
        :return: The sequence as a string.
        """
        return self.sequence

    def get_annotations(self) -> List[Annotation]:
        """
        Retrieves the complete list of Annotation objects currently associated with this part.
        :return: A list of `Annotation` objects. Returns an empty list if no annotations have been added.
        """
        return self.annotations

    def add_annotation(self, annotation: Annotation) -> AnnotatedPart:
        """
        Adds a single Annotation object to the current list of annotations for this part.
        :param annotation: The `Annotation` object to be added.
        :return: The `AnnotatedPart` instance itself, allowing for method chaining.
        """
        if self.annotations is None: # This check might be redundant if initialized to [] in __init__
            self.annotations = [annotation]
            return self

        self.annotations.append(annotation)

        return self

    def add_annotations(self, seq_annotations: List[Annotation]) -> AnnotatedPart:
        """
        Adds a sequence (list) of Annotation objects to the current list of annotations for this part.
        :param seq_annotations: A list of `Annotation` objects to be added.
        :return: The `AnnotatedPart` instance itself, allowing for method chaining.
        """
        if self.annotations is None: # This check might be redundant if initialized to [] in __init__
            self.annotations = seq_annotations
            return self

        self.annotations += seq_annotations

        return self

    def add(self, sequence: str, annotation: Optional[Annotation] = None) -> AnnotatedPart:
        """
        Appends a new sequence segment to the existing sequence of this part and
        optionally adds a new annotation.
        :param sequence: The sequence string to append.
        :param annotation: An optional `Annotation` object to add, corresponding to the appended sequence or a sub-region thereof.
        :return: The `AnnotatedPart` instance itself, allowing for method chaining.
        """
        if not self.sequence:
            self.sequence = ""

        if not self.annotations:
            self.annotations = []

        self.sequence += sequence
        if annotation:
            self.annotations.append(annotation)

        return self

    def get_seq_record(self) -> SeqRecord:
        """
        Converts the AnnotatedPart into a Biopython SeqRecord object.
        This object encapsulates the sequence, ID, name, description, and all associated features (annotations)
        in a format suitable for bioinformatics operations and file output.
        :return: A Biopython `Bio.SeqRecord.SeqRecord` object.
        """
        return SeqRecord(Seq(self.sequence), id=self.part_id, name=self.name, description=self.description, features=[
            annotation.to_seq_feature() for annotation in self.annotations], annotations={"molecule_type": "PROTEIN"})

    def to_genbank_string(self) -> str:
        """
        Generates a GenBank formatted string representation of this AnnotatedPart.
        :return: A string containing the GenBank record.
        """
        f: StringIO = StringIO()
        SeqIO.write(self.get_seq_record(), f, "genbank")
        return f.getvalue()

    def save_genbank_file(self, file_path: str) -> None:
        """
        Saves the GenBank formatted string of this AnnotatedPart to a specified file.
        :param file_path: The path to the file where the GenBank record will be saved.
        :return: None
        """
        content: str = self.to_genbank_string()
        with open(file_path, "w") as f:
            f.write(content)
