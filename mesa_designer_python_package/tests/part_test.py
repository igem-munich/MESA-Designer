from mesa_designer.mesa import *

def test_annotation():
    a = Annotation("test", 5, 10, part_type="CDS")

    ap = AnnotatedPart("ABCDEFGHIJKLMNOPQRSTUVWXYZ", name="test_part")

    ap.add_annotation(a)

    ap.save_genbank_file("test.gb")