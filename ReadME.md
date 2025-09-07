# iGEM Munich 2025 Mesa Design Tool

## Design Outline
[Design Doc](https://docs.google.com/document/d/1ciPsgLo5JNp7wKqFREEnWSCiShK2ZrRBRCNdM3VBm_A/edit?tab=t.0)

## Installation
### Install it locally
Start by installing python 3.13.x. You can get it from: https://www.python.org/downloads/ Check that you have the correct version installed by running: `python --version`

Next clone the github repo using any method you want, but we recommend running this command: `git clone https://github.com/igem-munich/MESA-Designer.git`

Next you will probably want to create a virtual environment. 
This can be done by running: `python -m venv .venv` This will be used to cleanly install all required packages. This can be done by running `pip install -r requirements.txt` from the cloned github repo's main directory.

As a final setup step you need to download all database files. This can be done by simply running the setup.py script: `python3 ./setup.py` from the main directory.

To start the tool, please run `streamlit run ./app/main.py` from the main directory!

### Build the docker image yourself

You can build the docker image with:
```bash
git clone https://gitlab.igem.org/2025/software-tools/munich mesa-designer
cd mesa-designer
docker build -t mesa-designer .
 ```
now tht the repository is build, run it
```bash
docker run -p 8501:8501 mesa-designer 
```
### Pulling image from dockerhub
The MESA-Designer tool can now also be run as a docker container. Please find it [here](https://hub.docker.com/repository/docker/aeneastews/mesa-designer/general).

#### Running docker container
You can run the docker container using : `docker run -p 8501:8501 mesa-designer`

It can then be accessed from a webbrowser at localhost:8501
