# iGEM Munich 2025 Mesa Design Tool

## Design Outline
[Design Doc](https://docs.google.com/document/d/1ciPsgLo5JNp7wKqFREEnWSCiShK2ZrRBRCNdM3VBm_A/edit?tab=t.0)

## Installation
Start by installing python 3.13.x. You can get it from: https://www.python.org/downloads/ Check that you have the correct version installed by running: `python --version`

Next clone the github repo using any method you want, but we recommend running this command: `git clone https://github.com/igem-munich/MESA-Designer.git`

Next you will probably want to create a virtual environment. 
This can be done by running: `python -m venv .venv` This will be used to cleanly install all required packages. This can be done by running `pip install -r requirements.txt` from the cloned github repo's main directory.

As a final setup step you need to download all necessary structures and database files. This can be done by simply running the setup.py script: `python3 ./setup.py` from the main directory.

## Usage
To start the tool, please run `streamlit run ./app/main.py` from the main directory!