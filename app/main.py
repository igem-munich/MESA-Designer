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
from util import TMD_DATA
from util.pdb_interaction import extract_fasta_from_pdb
import zipfile

# session state
if "sabdab" not in st.session_state:
    st.session_state.sabdab = None
if "skempi" not in st.session_state:
    st.session_state.skempi = None
if "pdbs" not in st.session_state:
    st.session_state.pdbs = None
if "tmd" not in st.session_state:
    st.session_state.tmd = None
if "linkers" not in st.session_state:
    st.session_state.linkers = {"linker1": None, "linker2": None}
if "pdb_fasta" not in st.session_state:
    st.session_state.pdb_fasta = None

st.title("MESA-Design Tool")

# columns for search field and search button
col1, col2 = st.columns([1, 0.1])

# create search field and button
with col1:
    search_field = st.text_input(label="Antigen-Search", key="search_input", label_visibility="collapsed", placeholder="Search target antigen")

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


# OPTIONAL: more options for linker building
if st.session_state.pdbs:
    st.header("Linker design")

    multi_linker = st.toggle(
        label="Different Linkers",
        value=False,
        key="multi_link",
        help="Build identical Linker for both chains or separate linkers"
    )

    # input area for 1 linker
    linker1_repeat_generation = st.container()
    linker1_columns = st.columns(3, vertical_alignment="bottom")
    with linker1_repeat_generation:
        with linker1_columns[0]:
            linker1_input = st.text_input(
                key="linker1_pattern",
                label="Linker Pattern",
                value="GGGS",
                max_chars=1000
            )
        with linker1_columns[1]:
            linker1_repeats = st.number_input(
                key="linker1_repeats",
                label="Repeats",
                min_value=1,
                max_value=1000,
                value=10
            )
        with linker1_columns[2]:
            linker1_pattern_generate = st.button(
                key="linker1_generate",
                label="Generate"
            )

            if linker1_pattern_generate:
                st.session_state.linkers["linker1"] = linker1_input * linker1_repeats

    def update_linker1():
        st.session_state.linkers["linker1"] = st.session_state.linker1_sequence

    linker1_sequence = st.text_input(
        key="linker1_sequence",
        label="Linker sequence",
        value=st.session_state.linkers["linker1"],
        on_change=update_linker1,
        max_chars=1000
    )

    # one linker design is always required, but second may be desired
    if multi_linker:
        linker2_repeat_generation = st.container()
        linker2_columns = st.columns(3, vertical_alignment="bottom")
        with linker2_repeat_generation:
            with linker2_columns[0]:
                linker2_input = st.text_input(
                    key="linker2_pattern",
                    label="Linker Pattern",
                    value="GGGS",
                    max_chars=1000
                )
            with linker2_columns[1]:
                linker2_repeats = st.number_input(
                    key="linker2_repeats",
                    label="Repeats",
                    min_value=1,
                    max_value=1000,
                    value=10
                )
            with linker2_columns[2]:
                linker2_pattern_generate = st.button(
                    key="linker2_generate",
                    label="Generate"
                )

                if linker2_pattern_generate:
                    st.session_state.linkers["linker2"] = linker2_input * linker2_repeats

        def update_linker2():
            st.session_state.linkers["linker2"] = st.session_state.linker2_sequence

        linker2_sequence = st.text_input(
            key="linker2_sequence",
            label="Linker sequence",
            value=st.session_state.linkers["linker2"],
            on_change=update_linker2,
            max_chars=1000
        )

# TMD picker
if st.session_state.pdbs:
    st.header("TMD Picker")

    multi_tmd = st.toggle(
        label="Different TMDs",
        value=False,
        key="multi_tmd",
        help="Pick identical TMD for both chains or separate TMDs"
    )

    if multi_tmd:
        col_tmd_left, col_tmd_right = st.columns(2)

        with col_tmd_left:
            st.subheader("TMD 1")
            selected_tmd_left = st.radio(
                "Select TMD 1",
                options=list(TMD_DATA.keys()),
                key="tmd_1",
                label_visibility="collapsed",
                index=0
            )

            st.code(TMD_DATA[selected_tmd_left][1], language="text")

        with col_tmd_right:
            st.subheader("TMD 2")
            selected_tmd_right = st.radio(
                "Select TMD 2",
                options=list(TMD_DATA.keys()),
                key="tmd_2",
                label_visibility="collapsed",
                index=4
            )

            st.code(TMD_DATA[selected_tmd_right][1], language="text")

        st.session_state.tmd = {"left": TMD_DATA[selected_tmd_left], "right": TMD_DATA[selected_tmd_right]}

    else:
        st.subheader("TMD 1,2")

        col_tmd_select, col_tmd_sequence = st.columns([0.5, 1], vertical_alignment="bottom")
        with col_tmd_select:
            selected_tmd = st.radio(
                "Selected TMD",
                options=list(TMD_DATA.keys()),
                key="tmd_both",
                label_visibility="collapsed",
                index=4
            )

        with col_tmd_sequence:
            st.write("Optional Display of TMD")
            st.code(TMD_DATA[selected_tmd][1], language="text")

        st.session_state.tmd = {"both": TMD_DATA[selected_tmd]}

# OPTIONAL: view the different TMDs and view their combinations
# TODO: intracellular
# TODO: complete download into zip file of all files
# download button for structure
#        with zipfile.ZipFile(file="download.zip")
#        with open(st.session_state.pdbs[pdb_selection], "rb") as f:
#            st.session_state.pdb_fasta = extract_fasta_from_pdb(st.session_state.pdbs[pdb_selection])
#            st.download_button(label="⬇️", data=f, file_name=f"{pdb_selection}.pdb")

print(st.session_state.linkers)