from mesa_designer import TMD_DATA, AIP_DATA, FRET_ICDs, CTEV_DATA, NTEV_DATA, TEVP_DATA, PRS_DATA, SIGNAL_SEQS, TAG_SEQS

def test_data_imports():
    assert len(TMD_DATA) > 0
    assert len(AIP_DATA) > 0
    assert len(FRET_ICDs) > 0
    assert len(CTEV_DATA) > 0
    assert len(NTEV_DATA) > 0
    assert len(TEVP_DATA) > 0
    assert len(PRS_DATA) > 0
    assert len(SIGNAL_SEQS) > 0
    assert len(TAG_SEQS) > 0
