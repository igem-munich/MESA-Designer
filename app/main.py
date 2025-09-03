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
from annotated_text import annotated_text
from util.antibody_search import search_antibodies
from util import TMD_DATA, CTEV_DATA, NTEV_DATA, TEVP_DATA, PRS_DATA, AIP_DATA, FRET_ICDs
from util.pdb_interaction import extract_chains_from_pdb
from util.general import new_random_color
import zipfile
import io
import json
from datetime import datetime

# session state
state = st.session_state

if "sabdab" not in state:
    state.sabdab = None
if "skempi" not in state:
    state.skempi = None
if "pdbs" not in state:
    state.pdbs = None
if "tmds" not in state:
    state.tmds = {}
if "linkers" not in state:
    state.linkers = {}
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
if "protease_chain_association" not in state:
    state.protease_chain_association = {"n": set(), "c": set(), "complete": set()}
if "target_chain_association" not in state:
    state.target_chain_association = set()
if "aip_chain_association" not in state:
    state.aip_chain_association = set()
if "construct_list_formatted" not in state:
    state.construct_list_formatted = {}
if "download_data" not in state:
    state.download_data = None

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

def update_linker_text_input(chain_id):
    state.linkers[f"{chain_id}_linker"] = state[f"{chain_id}_linker_sequence"].upper()

def generate_download() -> None:
    zip_buffer = io.BytesIO()

    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
       # add constructs
        for key, construct in state.construct_list_formatted.items():
            # skip unselected
            if not state[f"{key}_checkbox"]:
                continue

            text = construct[0].strip() + "\n"
            for part in construct[1:]:
                if isinstance(part, tuple):
                    text += part[0]
                else:
                    text += part

            file_name = f"{key.replace(' ', '_').replace(':', '').replace('-', '_').replace('/', '_')}.fasta"
            zf.writestr(file_name, text)

        # add PDB structures
        if state["download_sel_pdb"] and state.pdbs and state.prev_pdb_selection:
            selected_pdb_path = state.pdbs.get(state.prev_pdb_selection)

            if selected_pdb_path and Path(selected_pdb_path).exists():
                pdb_file_name = Path(selected_pdb_path).name

                try:
                    with open(selected_pdb_path, "rb") as f:
                        zf.writestr(f"selected_pdb/{pdb_file_name}", f.read())  # Put in subfolder
                except Exception as e:
                    st.warning(f"Could not add selected PDB '{pdb_file_name}' to zip: {e}")

        if state["download_all_pdb"] and state.pdbs:
            for pdb_id, pdb_path in state.pdbs.items():
                if pdb_path and Path(pdb_path).exists():
                    pdb_file_name = Path(pdb_path).name

                    try:
                        with open(pdb_path, "rb") as f:
                            zf.writestr(f"all_pdbs/{pdb_file_name}", f.read())  # Put in subfolder
                    except Exception as e:
                        st.warning(f"Could not add PDB '{pdb_file_name}' to zip: {e}")

        # add additional data
        if state["download_additional"]:  # Use the new checkbox key
            summary_content = f"""
                MESA Design Tool Output Summary
                ------------------------------
                Generated on: {datetime.today().strftime('%Y-%m-%d %H:%M:%S')}
                Selected PDB: {state.prev_pdb_selection if state.prev_pdb_selection else 'N/A'}

                Selected Binder FASTA:
                {state.pdb_fasta if state.pdb_fasta else 'N/A'}

                Linker Information:
                {json.dumps(state.linkers, indent=2)}

                TMD Information:
                {json.dumps(state.tmds, indent=2)}
                """
            zf.writestr("mesa_design_summary.txt", summary_content.strip())

            s = io.StringIO()
            state.sabdab.to_csv(s)
            zf.writestr("sabdab_data.csv", s.getvalue())

    zip_buffer.seek(0)
    state.download_data = zip_buffer.getvalue()


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
        state.pdb_fasta = fasta

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
if state.pdbs and len(state.highlight_selection) > 0:
    st.divider()
    st.header("Linker Design")

    if len(state.highlight_selection) < 1:
        st.error("Select chains to create linkers")

    # create linker creation menu for each selected chain
    for chain_id in state.highlight_selection:
        # update linker state
        if f"{chain_id}_linker" not in state.linkers:
            state.linkers[f"{chain_id}_linker"] = ""

        st.subheader(f"{chain_id} Linker")

        linker_columns = st.columns(3, vertical_alignment="bottom")
        with linker_columns[0]:
            linker_input = st.text_input(
                key=f"{chain_id}_linker_pattern",
                label=f"{chain_id} Linker Pattern",
                value="GGGS",
                max_chars=1000
            )
        with linker_columns[1]:
            linker_repeats = st.number_input(
                key=f"{chain_id}_linker_repeats",
                label="Repeats",
                min_value=1,
                max_value=1000,
                value=10
            )
        with linker_columns[2]:
            linker_pattern_generate = st.button(
                key=f"{chain_id}_linker_generate",
                label="Generate"
            )

            if state[f"{chain_id}_linker_generate"]:
                state.linkers[f"{chain_id}_linker"] = state[f"{chain_id}_linker_pattern"] * state[f"{chain_id}_linker_repeats"]

        linker_sequence = st.text_input(
            key=f"{chain_id}_linker_sequence",
            label=f"{chain_id} Linker Sequence",
            value=state.linkers[f"{chain_id}_linker"],
            on_change=update_linker_text_input,
            args=chain_id,
            max_chars=1000
        )

### TMD picker #########################################################################################################
# TODO: let user enter custom TMDs
if state.pdbs and len(state.highlight_selection) > 0:
    st.divider()
    st.header("TMD Picker")

    if len(state.highlight_selection) < 1:
        st.error("Select chains to create TMDS")

    for chain_id in state.highlight_selection:
        # update tmd state
        if f"{chain_id}_tmd" not in state.tmds:
            state.tmds[f"{chain_id}_tmd"] = ""

        st.subheader(f"{chain_id} TMD")

        tmd_selection = st.radio(
            label="tmd_selection",
            options=TMD_DATA.keys(),
            horizontal=True,
            label_visibility="collapsed",
            key=f"{chain_id}_tmd_selection"
        )
        st.code(
            TMD_DATA[tmd_selection][1],
            language="text",
            height="content",
            wrap_lines=True
        )

        state.tmds[f"{chain_id}_tmd"] = TMD_DATA[tmd_selection][1]


# TODO OPTIONAL: view the different TMDs and view their combinations
### INTRACELLULAR PART DESIGNER ########################################################################################
if state.pdbs and len(state.highlight_selection) > 0:
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
    if custom_icd:
        st.subheader("Enter ICD Sequence")
        custom_icd_sequence = st.text_area(
            label="Custom ICD",
            value="",
            height="content",
            label_visibility="collapsed",
            key="custom_icd_sequence"
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
                            label_visibility="collapsed",
                            key="n_protease_selection"
                        )

                        n_protease_sequence = st.code(
                            NTEV_DATA[n_protease_selection][1],
                            language="text",
                            wrap_lines=True,
                            height="content"
                        )

                        st.markdown("##### Append to Chain")
                        protease_association_cols_n = st.columns(len(state.highlight_selection))
                        for col, chain_id in zip(protease_association_cols_n, state.highlight_selection):
                            with col:
                                checkbox = st.checkbox(
                                    chain_id,
                                    key=f"{chain_id}_protease_association_n"
                                )

                                if checkbox:
                                    state.protease_chain_association["n"].add(chain_id)
                                else:
                                    if chain_id in state.protease_chain_association["n"]:
                                        state.protease_chain_association["n"].remove(chain_id)

                    with split_c_col:
                        st.markdown("#### C-Terminus")

                        # split protease selection
                        c_protease_selection = st.radio(
                            "Pick C-Terminal TEV-Protease",
                            options=CTEV_DATA.keys(),
                            label_visibility="collapsed",
                            key="c_protease_selection"
                        )

                        c_protease_sequence = st.code(
                            CTEV_DATA[c_protease_selection][1],
                            language="text",
                            wrap_lines=True,
                            height="content"
                        )

                        st.markdown("##### Append to Chain")
                        protease_association_cols_c = st.columns(len(state.highlight_selection))
                        for col, chain_id in zip(protease_association_cols_c, state.highlight_selection):
                            with col:
                                checkbox = st.checkbox(
                                    chain_id,
                                    key=f"{chain_id}_protease_association_c"
                                )

                                if checkbox:
                                    state.protease_chain_association["c"].add(chain_id)
                                else:
                                    if chain_id in state.protease_chain_association:
                                        state.protease_chain_association["c"].remove(chain_id)
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
                            label_visibility="collapsed",
                            key="n_protease_sequence_entry"
                        )

                    with split_c_col:
                        st.markdown("##### C-Terminus")
                        c_protease_sequence_entry = st.text_area(
                            label="C-Terminus",
                            value="",
                            height="content",
                            max_chars=5000,
                            placeholder="Enter c-terminal protease sequence",
                            label_visibility="collapsed",
                            key="c_protease_sequence_entry"
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
                            label_visibility="collapsed",
                            key="complete_protease_selection"
                        )

                    with sequence_col:
                        st.code(
                            TEVP_DATA[complete_protease_selection][1],
                            language="text",
                            wrap_lines=True,
                            height="content"
                        )

                    st.markdown("##### Append to Chain")
                    protease_association_cols = st.columns(len(state.highlight_selection))
                    for col, chain_id in zip(protease_association_cols, state.highlight_selection):
                        with col:
                            checkbox = st.checkbox(
                                chain_id,
                                key=f"{chain_id}_protease_association"
                            )

                            if checkbox:
                                state.protease_chain_association["complete"].add(chain_id)
                            else:
                                if chain_id in state.protease_chain_association:
                                    state.protease_chain_association["complete"].remove(chain_id)
                else:
                    st.markdown("#### Protease Sequence")

                    complete_protease_sequence = st.text_area(
                        label="Protease Sequence",
                        value="",
                        height="content",
                        max_chars=5000,
                        label_visibility="collapsed",
                        placeholder="Enter protease sequence",
                        key="custom_protease_sequence"
                    )

                # let user pick which binder chains the protease should be attached to


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
                        label_visibility="collapsed",
                        key="prs_selection"
                    )

                with prs_seq_col:
                    prs_sequence_display = st.code(
                        PRS_DATA[prs_selection][1],
                        language="text",
                        wrap_lines=True,
                        height="content"
                    )

            # enter custom PRS sequence
            else:
                st.markdown("#### PRS Sequence")
                prs_sequence = st.text_area(
                    label="Custom PRS Sequence",
                    value="",
                    height="content",
                    max_chars=5000,
                    label_visibility="collapsed",
                    key="custom_prs_sequence",
                    placeholder="Enter PRS sequence"
                )

            # target sequence
            st.markdown("#### Target Sequence")
            target_sequence = st.text_area(
                label="Target Sequence",
                value="",
                height="content",
                label_visibility="collapsed",
                placeholder="Enter target sequence",
                max_chars=10000,
                key="target_sequence"
            )

            # associate target to specific chains
            st.markdown("##### Append to Chain")
            target_chain_cols = st.columns(len(state.highlight_selection))
            for col, chain_id in zip(target_chain_cols, state.highlight_selection):
                with col:
                    checkbox = st.checkbox(
                        label=chain_id,
                        key=f"{chain_id}_target_association"
                    )

                    if checkbox:
                        state.target_chain_association.add(chain_id)
                    else:
                        if chain_id in state.target_chain_association:
                            state.target_chain_association.remove(chain_id)

    # further options
    st.subheader("Further Options")
    further_options_container = st.container(border=True)

    with further_options_container:
        if not custom_icd:
            # include auto-inhibitory peptide
            aip_toggle = st.toggle(
                label="Include AIP",
                value=False,
                key="include_aip",
                help="Include Auto-Inhibitory Peptide. Can reduce background noise by inhibiting Protease"
            )

            if aip_toggle:
                custom_aip = st.toggle(
                    label="Custom AIP",
                    value=False,
                    key="custom_aip_toggle",
                    help="Choose from variants of TEVp AIPs or create custom"
                )

                # pick TEVp aip
                if not custom_aip:
                    st.markdown("#### AIP Selection")
                    aip_selection_col, aip_seq_col = st.columns([0.5, 2], vertical_alignment="bottom")
                    with aip_selection_col:
                        aip_selection = st.radio(
                            label="AIP Selection",
                            options=AIP_DATA.keys(),
                            label_visibility="collapsed",
                            key="aip_selection"
                        )

                    with aip_seq_col:
                        st.code(
                            AIP_DATA[aip_selection][1],
                            language="text",
                            wrap_lines=True,
                            height="content"
                        )

                # custom AIP
                else:
                    st.markdown("#### Custom AIP")

                    aip_sequence = st.text_area(
                        label="AIP Sequence",
                        value="",
                        height="content",
                        max_chars=1000,
                        label_visibility="collapsed",
                        placeholder="Enter AIP sequence",
                        key="custom_aip_sequence"
                    )

                # associate AIP with chains
                st.markdown("##### Append to Chain")
                aip_chain_cols = st.columns(len(state.highlight_selection))
                for col, chain_id in zip(aip_chain_cols, state.highlight_selection):
                    with col:
                        checkbox = st.checkbox(
                            chain_id,
                            key=f"{chain_id}_aip_association"
                        )

                        if checkbox:
                            state.aip_chain_association.add(chain_id)
                        else:
                            if chain_id in state.aip_chain_association:
                                state.aip_chain_association.remove(chain_id)

            else:
                state.aip_chain_association.clear()

            # options divider
            st.divider()

        # automatically create FRET chains
        fret_chains_toggle = st.toggle(
            label="Create FRET sequences",
            value=False,
            key="fret_chains_toggle",
            help="Create FRET chains to test selected binders and TMDs"
        )

### DOWNLOAD AND OVERVIEW ##############################################################################################
if state.pdbs and len(state.highlight_selection) > 0:
    st.divider()
    st.header("Component Overview")
    overview_container = st.container(border=True)

    with overview_container:
        # sections
        # binder fasta display
        st.subheader("Binder Overview")
        st.code(
            state.pdb_fasta,
            language="text",
            wrap_lines=True,
            height="content"
        )

        # linker combination
        st.divider()
        st.subheader(f"ECD Linker{'s' if len(state.highlight_selection) > 1 else ''}")

        for chain_id in state.highlight_selection:
            st.markdown(f"#### {chain_id} Linker")
            st.code(
                state.linkers[f"{chain_id}_linker"],
                language="text",
                wrap_lines=True,
                height="content"
            )

        # tmd overview
        st.divider()
        st.subheader("TMD Overview")

        for chain_id in state.highlight_selection:
            st.markdown(f"#### {chain_id} TMD")
            st.code(
                state.tmds[f"{chain_id}_tmd"],
                language="text",
                wrap_lines=True,
                height="content"
            )


        # icd overview
        st.divider()
        if state.custom_icd_toggle:
            st.subheader("Custom ICD")
            st.code(
                state.custom_icd_sequence,
                language="text",
                wrap_lines=True,
                height="content"
            )
        else:
            if state.split_protease_toggle:
                st.subheader("Split Protease Overview")
                split_protease_overview_cols = st.columns(2)

                with split_protease_overview_cols[0]:
                    st.markdown("#### N-Terminus")
                    st.code(
                        NTEV_DATA[state.n_protease_selection][1] if not state.custom_protease_toggle else state.n_protease_sequence_entry,
                        language="text",
                        wrap_lines=True,
                        height="content"
                    )

                with split_protease_overview_cols[1]:
                    st.markdown("#### C-Terminus")
                    st.code(
                        CTEV_DATA[state.c_protease_selection][1] if not state.custom_protease_toggle else state.c_protease_sequence_entry,
                        language="text",
                        wrap_lines=True,
                        height="content"
                    )

            else:
                st.subheader("Protease Overview")
                st.markdown("#### Protease Sequence")
                st.code(
                    TEVP_DATA[state.complete_protease_selection][1] if not state.custom_protease_toggle else state.custom_protease_sequence,
                    language="text",
                    wrap_lines=True,
                    height="content"
                )

            if state.include_aip:
                st.markdown(f"#### {'Custom ' if state.custom_aip_toggle else ' '}AIP Sequence")
                st.code(
                    AIP_DATA[state.aip_selection][1] if not state.custom_aip_toggle else state.custom_aip_sequence,
                    language="text",
                    wrap_lines=True,
                    height="content"
                )


            st.divider()
            st.subheader("Target Overview")
            st.markdown(f"#### {'Custom ' if state.custom_prs_toggle else ' '}PRS Sequence")
            st.code(
                PRS_DATA[state.prs_selection][1] if not state.custom_prs_toggle else state.custom_prs_sequence,
                wrap_lines=True,
                language="text",
                height="content"
            )

            st.markdown("#### Target Sequence")
            st.code(
                state.target_sequence,
                language="text",
                wrap_lines=True,
                height="content"
            )

    st.header("Downloads")
    download_container = st.container(border=True)
    with download_container:
        # assemble annotated constructs
        construct_list: dict[str, list[str | tuple[str, str] | tuple[str, str, str]]] = {}
        for chain_id in state.highlight_selection:
            for chain_data in state.current_pdb_chains_data:
                if chain_data["chain_id"] != chain_id:
                    continue

                current_chain: list[str | tuple[str, str] | tuple[str, str, str]] = [
                    (chain_data["sequence"], "Binder", "#534cb3"),
                    (state.linkers[f"{chain_id}_linker"], "Linker", "#eba814"),
                    (state.tmds[f"{chain_id}_tmd"], "TMD", "#69ad52"),
                ]

                # selectively append target and protease
                if chain_id in state.target_chain_association:
                    construct_list[f"{chain_id}_Target"] = ([f"> {chain_id}_Target\n\n"] + current_chain + [(PRS_DATA[state.prs_selection][1] if not state.custom_prs_toggle else state.custom_prs_sequence, "PRS", "#b4774b"), (state.target_sequence, "Target", "#bd4258")])

                if state.split_protease_toggle:
                    if chain_id in state.protease_chain_association["n"]:
                        construct_list[f"{chain_id}_N-Term Protease"] = ([f"> {chain_id}_N-Term Protease\n\n"] + current_chain + [(NTEV_DATA[state.n_protease_selection][1] if not state.custom_protease_toggle else state.n_protease_sequence_entry, "N-Term Protease", "#bfbd40")])
                        if chain_id in state.aip_chain_association:
                            construct_list[f"{chain_id}_N-Term Protease"].append((AIP_DATA[state.aip_selection][1] if not state.custom_aip_toggle else state.custom_aip_sequence, "AIP", "#5aa56b"))
                    elif chain_id in state.protease_chain_association["c"]:
                        construct_list[f"{chain_id}_C-Term Protease"] = ([f"> {chain_id}_C-Term Protease\n\n"] + current_chain + [(CTEV_DATA[state.c_protease_selection][1] if not state.custom_protease_toggle else state.c_protease_sequence_entry, "C-Term Protease", "#3948c6")])
                        if chain_id in state.aip_chain_association:
                            construct_list[f"{chain_id}_C-Term Protease"].append((AIP_DATA[state.aip_selection][1] if not state.custom_aip_toggle else state.custom_aip_sequence, "AIP", "#5aa56b"))

                else:
                    if chain_id in state.protease_chain_association["complete"]:
                        construct_list[f"{chain_id}_Protease"] = ([f"> {chain_id}_Protease\n\n"] + current_chain + [(TEVP_DATA[state.complete_protease_selection][1] if not state.custom_protease_toggle else state.custom_protease_sequence, "Protease", "#6b46b9")])
                        if chain_id in state.aip_chain_association:
                            construct_list[f"{chain_id}_Protease"].append((AIP_DATA[state.aip_selection][1] if not state.custom_aip_toggle else state.custom_aip_sequence, "AIP", "#5aa56b"))

                # create FRET sequences
                if state.fret_chains_toggle:
                    construct_list[f"{chain_id}_FRET_mVenus"] = ([f"> {chain_id}_FRET_mVenus\n\n"] + current_chain + [(FRET_ICDs["mVenus"][1], "mVenus", "#43b6bc")])
                    construct_list[f"{chain_id}_FRET_mCerulean"] = ([f"> {chain_id}_FRET_mCerulean\n\n"] + current_chain + [(FRET_ICDs["mCerulean"][1], "mCerulean", "#c43b81")])

        # store in state
        state.construct_list_formatted = construct_list

        # display constructs
        # construct and download selection
        st.subheader("Construct Overview")

        for key, construct in construct_list.items():
            with st.container():
                download_cols = st.columns([4, 0.1], vertical_alignment="center")
                with download_cols[0]:
                    with st.container(border=True):
                        annotated_text(construct)

                with download_cols[1]:
                    checkbox = st.checkbox(
                        label="Select for download",
                        value=True,
                        key=f"{key}_checkbox",
                        label_visibility="collapsed"
                    )

        # include additional files such as pdb and data
        st.subheader("Additional Files")
        with st.container(border=True):
            additional_cols = st.columns(2)
            with additional_cols[0]:
                st.checkbox(
                    label="Include PDB Structure",
                    value=True,
                    key="download_sel_pdb"
                )

                st.checkbox(
                    label="Include all PDB Structures",
                    value=False,
                    key="download_all_pdb"
                )

            with additional_cols[1]:
                st.checkbox(
                    label="Include additional Data",
                    value=True,
                    key="download_additional"
                )

        # download button
        if state.construct_list_formatted:
            generate_download()

        if state.download_data:
            download = st.download_button(
                label="Download Selected",
                key="download_selected",
                icon=":material/download:",
                data=state.download_data,
                file_name="mesa-design.zip",
                mime="application/zip"
            )

# TODO: unify naming conventions
# TODO: provide comments
# TODO: fasta annotation from uniprot
# TODO: genbank downloads
# TODO: score every component
# TODO: LLM questions
# TODO: docker container