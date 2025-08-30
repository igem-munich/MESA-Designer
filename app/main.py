from pathlib import Path
import sys
# add main directory to system path for imports
current_dir = Path(__file__).resolve().parent
project_root = current_dir.parent
sys.path.insert(0, str(project_root))

# regular code
import streamlit as st
from stmol import *
import py3Dmol
from util.antibody_search import search_antibodies
from util import TMD_DATA, CTEV_DATA, NTEV_DATA, TEVP_DATA, PRS_DATA, AIP_DATA
from util.pdb_interaction import extract_chains_from_pdb
from util.general import new_random_color
import zipfile

# session state
state = st.session_state

if "sabdab" not in state:
    state.sabdab = None
if "skempi" not in state:
    state.skempi = None
if "pdbs" not in state:
    state.pdbs = None
if "tmd" not in state:
    state.tmd = None
if "linkers" not in state:
    state.linkers = {"linker1": " ", "linker2": " "}
if "pdb_fasta" not in state:
    state.pdb_fasta = None
if "current_pdb_chains_data" not in state:
    state.current_pdb_chains_data = []
if "highlight_selection" not in state:
    state.highlight_selection = {}
if "prev_pdb_selection" not in state:
    state.prev_pdb_selection = None
if "chain_colors" not in state:
    state.chain_colors = {}
if "binder_fasta" not in state:
    state.binder_fasta = ""

# Callback function to update highlight_selection when a checkbox changes
def update_chain_highlight_selection(chain_id_to_toggle, current_pdb_selection):
    # Access the actual state of the specific checkbox that was changed
    checkbox_key = f"{current_pdb_selection}_checkbox_chain_{chain_id_to_toggle}"
    pdb_chain_data = extract_chains_from_pdb(state.pdbs[pdb_selection])
    lengths = {chain_data["chain_id"]: len(chain_data["sequence"]) for chain_data in pdb_chain_data}

    if state[checkbox_key]:  # If the checkbox is now checked
        if chain_id_to_toggle not in state.highlight_selection:
            state.highlight_selection[chain_id_to_toggle] = list(range(1, lengths[chain_id_to_toggle] + 1))
    else:  # If the checkbox is now unchecked
        if chain_id_to_toggle in state.highlight_selection:
            del state.highlight_selection[chain_id_to_toggle]

def update_chain_highlight_selection_residues(chain_id_to_change, current_pdb_selection):
    input_key = f"{current_pdb_selection}_residue_input_chain_{chain_id_to_change}"

    if chain_id_to_change in state.highlight_selection.keys():
        select_from, select_to = state[input_key].split(":")
        state.highlight_selection[chain_id_to_change] = list(range(int(select_from), int(select_to) + 1))

def update_split_protease_value():
    state.split_protease_toggle_value = not state.split_protease_toggle_value

# streamlit config
st.title("MESA-Design Tool")
st.set_page_config(layout="wide")

# columns for search field and search button
col1, col2 = st.columns([1, 0.1])

### Target Search Field ################################################################################################
# create search field and button
with col1:
    search_field = st.text_input(label="Antigen-Search", key="search_input", label_visibility="collapsed", placeholder="Search target antigen")

with col2:
    search_button = st.button("", key="search_button", icon=":material/search:", width=45)

# search database and display options
if search_button:
    if search_field:
        with st.spinner(f"Searching for: **{search_field}**"):
            state.sabdab, state.skempi, state.pdbs = search_antibodies(search_field)
    else:
        st.error("Please enter a search query.")

if state.sabdab is not None:
    if len(state.sabdab) > 0:
        st.dataframe(state.sabdab)
    else:
        st.info("No targets were found")

### Display Found Binder Structures ####################################################################################
# create new columns for displaying pdb structure
if state.pdbs:
    st.header("Inspect Structures")
    col3, col4 = st.columns([0.2, 1])

    # display radio selection and chain/residue selection
    with col3:

        st.subheader("Select Binder")
        pdb_selection = st.radio(
            label="PDB Selection",
            options=state.pdbs,
            label_visibility="collapsed"
        )

        # Re-extract chains and reset selection if PDB selection changes
        if state.prev_pdb_selection != pdb_selection:
            state.current_pdb_chains_data = extract_chains_from_pdb(state.pdbs[pdb_selection])
            state.highlight_selection = {}  # Clear previous selection

            # generate color for every chain
            state.chain_colors = {}
            for chain_data in state.current_pdb_chains_data:
                state.chain_colors[chain_data["chain_id"]] = new_random_color(list(state.chain_colors.values()))

            state.prev_pdb_selection = pdb_selection

        # Chain selection
        st.subheader("Select Chains")

        # Create checkboxes for each chain
        for chain_info in state.current_pdb_chains_data:
            chain_id = chain_info["chain_id"]

            # create columns for selection of chains and subsequent residue selection
            col_pdb_selection_left, col_pdb_selection_right = st.columns(2, vertical_alignment="bottom")
            with col_pdb_selection_left:
                st.checkbox(
                    f"Chain {chain_id}",
                    key=f"{pdb_selection}_checkbox_chain_{chain_id}",
                    value=(chain_id in state.highlight_selection.keys()),
                    on_change=update_chain_highlight_selection,
                    args=(chain_id, pdb_selection)
                )

            with col_pdb_selection_right:
                checkbox_key = f"{pdb_selection}_checkbox_chain_{chain_id}"
                if state[checkbox_key]:
                    st.text_input(
                        f"Residues in Chain {chain_id}",
                        key=f"{pdb_selection}_residue_input_chain_{chain_id}",
                        value=f"{state.highlight_selection[chain_id][0]}:{state.highlight_selection[chain_id][-1]}",
                        on_change=update_chain_highlight_selection_residues,
                        args=(chain_id, pdb_selection),
                        label_visibility="collapsed"
                    )

    with col4:
        # display selected residues using py3dmol
        with open(state.pdbs[pdb_selection], "r") as f:
            # create py3Dmol view
            view = py3Dmol.view(width=1200, height=500)

            # add protein model in cartoon view
            add_model(view, xyz=f.read(), model_style="cartoon")

        # set viewstyle based on selection
        for chain_id in state.highlight_selection.keys():
            for residue_index in state.highlight_selection[chain_id]:
                view.setStyle({"chain": chain_id, "resi": residue_index}, {"cartoon": {"color": state.chain_colors[chain_id], "arrows": True}})

        # TODO: extract non-protein sequences and add them separately in stick view (issue with finding small molecule)

        # add hover functionality (chain, residue, residue number)
        view.setHoverable({}, True,
                          '''function(atom,viewer,event,container) {
                if(!atom.label) {
                atom.label = viewer.addLabel(`Chain ${atom.chain}:${atom.resn}:${atom.resi}`,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
               }}''',
           '''function(atom,viewer) { 
               if(atom.label) {
                viewer.removeLabel(atom.label);
                delete atom.label;
               }
            }''')
        view.zoomTo()
        showmol(view, height=500, width=1200)

    # show current fasta selection to user
    if len(state.highlight_selection) != 0:
        fasta = ""

        # get data
        pdb_data = extract_chains_from_pdb(state.pdbs[pdb_selection])

        # iterate through pdb entries and add selected chains/residues to fasta string
        for chain_data in pdb_data:
            if chain_data["chain_id"] in state.highlight_selection.keys():
                chain_id = chain_data["chain_id"]
                selection = state.highlight_selection[chain_id]
                fasta += f">{chain_id}\n"
                fasta += f"{chain_data['sequence'][selection[0] - 1:selection[-1]]}\n"

        # remove last space
        fasta = fasta[:-1]

        # show fasta
        st.subheader("FASTA Preview")
        fasta_preview = st.code(
            body=fasta,
            height="content",
            language="text",
            wrap_lines=True
        )

        state.binder_fasta = fasta

### Build linker between Binder and TMD ################################################################################
# TODO OPTIONAL: more options for linker building
if state.pdbs:
    st.divider()
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
                state.linkers["linker1"] = state.linker1_pattern * state.linker1_repeats

    def update_linker1():
        state.linkers["linker1"] = state.linker1_sequence

    linker1_sequence = st.text_input(
        key="linker1_sequence",
        label="Linker sequence",
        value=state.linkers["linker1"],
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
                    state.linkers["linker2"] = linker2_input * linker2_repeats

        def update_linker2():
            state.linkers["linker2"] = state.linker2_sequence

        linker2_sequence = st.text_input(
            key="linker2_sequence",
            label="Linker sequence",
            value=state.linkers["linker2"],
            on_change=update_linker2,
            max_chars=1000
        )

### TMD picker #########################################################################################################
# TODO: let user enter custom TMDs
if state.pdbs:
    st.divider()
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

            st.code(TMD_DATA[selected_tmd_left][1], language="text", wrap_lines=True)

        with col_tmd_right:
            st.subheader("TMD 2")
            selected_tmd_right = st.radio(
                "Select TMD 2",
                options=list(TMD_DATA.keys()),
                key="tmd_2",
                label_visibility="collapsed",
                index=4
            )

            st.code(TMD_DATA[selected_tmd_right][1], language="text", wrap_lines=True)

        state.tmd = {"left": TMD_DATA[selected_tmd_left], "right": TMD_DATA[selected_tmd_right]}

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
            st.code(TMD_DATA[selected_tmd][1], language="text", wrap_lines=True)

        state.tmd = {"both": TMD_DATA[selected_tmd]}

# TODO OPTIONAL: view the different TMDs and view their combinations
# TODO: intracellular: receptor test: FRET
### INTRACELLULAR PART DESIGNER ########################################################################################
if state.pdbs:
    st.divider()
    st.header("Intracellular Component")

    # choose between custom and tev protease ICD
    custom_icd = st.toggle(
        label="Custom ICD",
        value=False,
        key="custom_icd_toggle",
        help="You can design your own custom ICD if you want or continue building with the TEV-Protease. Note: If you would like to use the general MESA ICD framework, but use a different protein / recognition sequence you can still continue to use this tool."
    )

    # custom icd design
    if not custom_icd:
        st.subheader("Enter ICD Sequence")
        custom_icd_sequence = st.text_area(
            label="Custom ICD",
            value="",
            height="content",
            label_visibility="collapsed"
        )

    # TEV-Protease ICD Design
    else:
        st.subheader("Protease Design")
        protease_design_container = st.container(border=True)

        with protease_design_container:

            # split protease design or separate chain design
            # add toggle value to state to change label
            if "split_protease_toggle_value" not in state:
                state.split_protease_toggle_value = True

            # pick between split protease or separate sequences
            split_protease = st.toggle(
                label=("Split Protease" if state.split_protease_toggle_value else "Separate Chains"),
                value=state.split_protease_toggle_value,
                key="split_protease_toggle",
                help="Split protease: The TEV-Protease or custom protein is split into two chains. Separate Chains: The TEV-Protease is fully assembled on one chain and the protein recognition sequence and target are on the other chain",
                on_change=update_split_protease_value
            )

            # let user pick between prebuilt TEV-protease construct or custom
            custom_protease = st.toggle(
                label="Custom Protease",
                value=False,
                key="custom_protease_toggle",
                help="Enter your own protease or use TEV-Protease"
            )

            # split protease design
            if state.split_protease_toggle_value:
                if not custom_protease:
                    split_n_col, split_c_col = st.columns(2)

                    # choose split parts from either pre-defined or custom
                    with split_n_col:
                        st.markdown("#### N-Terminus")

                        # split protease selection
                        n_protease_selection = st.radio(
                            "Pick N-Terminal TEV-Protease",
                            options=NTEV_DATA.keys(),
                            label_visibility="collapsed"
                        )

                        n_protease_sequence = st.code(
                            NTEV_DATA[n_protease_selection][1],
                            language="text",
                            wrap_lines=True
                        )

                    with split_c_col:
                        st.markdown("#### C-Terminus")

                        # split protease selection
                        c_protease_selection = st.radio(
                            "Pick C-Terminal TEV-Protease",
                            options=CTEV_DATA.keys(),
                            label_visibility="collapsed"
                        )

                        c_protease_sequence = st.code(
                            CTEV_DATA[c_protease_selection][1],
                            language="text",
                            wrap_lines=True
                        )
                else:
                    st.markdown("#### Custom Protease")
                    st.info("You can use [SPELL](https://dokhlab.med.psu.edu/spell/login.php) to guide your splitting process!")

                    split_n_col, split_c_col = st.columns(2)

                    with split_n_col:
                        st.markdown("##### N-Terminus")
                        n_protease_sequence_entry = st.text_area(
                            label="N-Terminus",
                            value="",
                            height="content",
                            max_chars=5000,
                            placeholder="Enter n-terminal protease sequence",
                            label_visibility="collapsed"
                        )

                    with split_c_col:
                        st.markdown("##### C-Terminus")
                        c_protease_sequence_entry = st.text_area(
                            label="C-Terminus",
                            value="",
                            height="content",
                            max_chars=5000,
                            placeholder="Enter c-terminal protease sequence",
                            label_visibility="collapsed"
                        )

            # complete protease design
            else:
                if not custom_protease:
                    st.markdown("#### Protease Selection")

                    selection_col, sequence_col = st.columns([0.5, 2], vertical_alignment="bottom")

                    with selection_col:
                        # select from existing proteases
                        complete_protease_selection = st.radio(
                            label="Select protease",
                            options=TEVP_DATA.keys(),
                            label_visibility="collapsed"
                        )

                    with sequence_col:
                        st.code(
                            TEVP_DATA[complete_protease_selection][1],
                            language="text",
                            wrap_lines=True
                        )
                else:
                    st.markdown("#### Protease Sequence")

                    complete_protease_sequence = st.text_area(
                        label="Protease Sequence",
                        value="",
                        height="content",
                        max_chars=5000,
                        label_visibility="collapsed",
                        placeholder="Enter protease sequence"
                    )

        # target design
        st.subheader("Target Design")
        target_design_container = st.container(border=True)

        with target_design_container:
            custom_prs = st.toggle(
                label="Custom PRS",
                value=False,
                key="custom_prs_toggle",
                help="Enter a custom PRS or use an existing version",
            )

            # pick from default TEVp PRS
            if not custom_prs:
                st.markdown("#### TEVp PRS Variant Selection")

                prs_selection_col, prs_seq_col = st.columns([0.5, 2], vertical_alignment="bottom")
                with prs_selection_col:
                    prs_selection = st.radio(
                        label="PRS Selection",
                        options=PRS_DATA.keys(),
                        label_visibility="collapsed"
                    )

                with prs_seq_col:
                    prs_sequence_display = st.code(
                        PRS_DATA[prs_selection][1],
                        language="text",
                        wrap_lines=True
                    )

            # enter custom PRS sequence
            else:
                st.markdown("#### PRS Sequence")
                prs_sequence = st.text_area(
                    label="Custom PRS Sequence",
                    value="",
                    height="content",
                    max_chars=5000,
                    label_visibility="collapsed"
                )

            # target sequence
            st.markdown("#### Target Sequence")
            target_sequence = st.text_area(
                label="Target Sequence",
                value="",
                height="content",
                label_visibility="collapsed",
                placeholder="Enter target sequence",
                max_chars=10000
            )

    # further options
    st.subheader("Further Options")
    further_options_container = st.container(border=True)

    with further_options_container:
        # include auto-inhibitory peptide
        aip_toggle = st.toggle(
            label="Include AIP",
            value=False,
            key="include_aip",
            help="Include Auto-Inhibitory Peptide. Can reduce background noise by inhibiting Protease"
        )

        if aip_toggle:
            st.markdown("#### AIP Selection")
            aip_selection_col, aip_seq_col = st.columns([0.5, 2], vertical_alignment="bottom")
            with aip_selection_col:
                aip_selection = st.radio(
                    label="AIP Selection",
                    options=AIP_DATA.keys(),
                    label_visibility="collapsed"
                )

            with aip_seq_col:
                st.code(
                    AIP_DATA[aip_selection][1],
                    language="text",
                    wrap_lines=True
                )

    # pick if ECD FRET test should be done after building complete construct




# TODO: TEV picker
#   - Protease/Target Chain or Split
#   - Protease/Target: pick protease and optional AIP && build target chain: pick PRS & enter target
#   - Split: link to SPELL for user to split their target protein or use default





# TODO: complete download into zip file of all files
# download button for structure
#        with zipfile.ZipFile(file="download.zip")
#        with open(state.pdbs[pdb_selection], "rb") as f:
#            state.pdb_fasta = extract_fasta_from_pdb(state.pdbs[pdb_selection])
#            st.download_button(label="⬇️", data=f, file_name=f"{pdb_selection}.pdb")

# TODO: use containers to create clear design units