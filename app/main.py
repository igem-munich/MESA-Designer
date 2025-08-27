from pathlib import Path
import sys
# add main directory to system path for imports
current_dir = Path(__file__).resolve().parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

# regular code
import streamlit as st
import streamlit_molstar as molstar
from util.antibody_search import search_antibodies

# session state
if "sabdab" not in st.session_state:
    st.session_state.sabdab = None
if "skempi" not in st.session_state:
    st.session_state.skempi = None
if "pdbs" not in st.session_state:
    st.session_state.pdbs = None

st.title("MESA-Design Tool")

# columns for search field and search button
col1, col2 = st.columns([1, 0.1])

# create search field and button
with col1:
    search_field = st.text_input(label="Antigen-Search", key="search_input", label_visibility="collapsed", placeholder="Enter your target Antigen")

with col2:
    search_button = st.button("🔎", key="search_button")

# search database and display options
if search_button:
    if search_field:
        st.write(f"Searching for: **{search_field}**")
        st.session_state.sabdab, st.session_state.skempi, st.session_state.pdbs = search_antibodies(search_field)
    else:
        st.error("Please enter a search query.")

if st.session_state.sabdab is not None:
    st.dataframe(st.session_state.sabdab)

# create new columns for displaying pdb structure
if st.session_state.pdbs:
    st.header("Inspect Structures")
    col3, col4 = st.columns([0.2, 1])

    # display radio selection
    with col3:
        pdb_selection = st.radio(
            label="PDB Selection",
            options=st.session_state.pdbs,
            label_visibility="collapsed"
        )

    with col4:
        molstar.st_molstar(st.session_state.pdbs[pdb_selection], height=500)


# TODO: create composition and length selection for linker
# TODO: TMD picker
# TODO: intracellular
