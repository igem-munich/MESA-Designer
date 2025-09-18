from mesa_designer.mesa import *

def test_mesa_chain():
    chain_a = MesaChain()
    chain_b = MesaChain(name="b")

    chain_a.add_binder(sequence="BINDERASEQUENCE")
    chain_b.add_binder(sequence="BINDERBSEQUENCE", name="b-binder", annotation="binder_b_annotation")

    chain_a.add_tmd_linker()
    chain_b.add_tmd_linker(sequence="BLINKER")

    chain_a.add_protease(protease_name=list(NTEV_DATA.keys())[0])
    chain_b.add_custom_protease(sequence="CUSTOMCPROTEASE", name="C Protease", annotation="c_protease_annotation")

    chain_a.add_prs(prs_name="PRS")
    chain_b.add_custom_prs(sequence="CUSTOMPRS", name="custom_b_prs", annotation="custom_b_annotation")

    chain_a.add_cargo(sequence="CHAINATARGET", name="chain_a_target", annotation="chain_a_target_annotation")

    chain_a.add_tmd("FGFR4")

    chain_a.add_signal_peptide(peptide_name="CD4")

    a_ap = chain_a.to_annotated_part(name="A_annotated_part")
    b_ap = chain_b.to_annotated_part(name="B_annotated_part")

    a_ap.save_genbank_file("a_ap.gb")
    b_ap.save_genbank_file("b_ap.gb")
    chain_a.save_genbank_file("a_chain.gb")
    chain_b.save_genbank_file("b_chain.gb")

    fret_assembly: MesaAssembly = chain_a.to_fret_chains()
    fret_assembly.save_genbank_files(".")