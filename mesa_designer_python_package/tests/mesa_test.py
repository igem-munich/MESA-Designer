from mesa_designer.mesa import *

chain_a = MesaChain()

chain_a.add_binder(sequence="BINDERASEQUENCE")

chain_a.save_genbank_file("a.gb")