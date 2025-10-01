# Team Munich 2025 Software Tool

## Description
The iGEM Munich 2025 MESA Designer Software Tool is an easy-to-use [MESA (Modular Extracellular Sensor Architecture for Engineering Mammalian Cell-based Devices)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4161666/) design tool.
The MESA Receptor Framework is a highly flexible and adaptable method for custom cellular circuit creation. The MESA Designer Tool provides a ... package for creation. Docker, Python package, webapp.
  
TODOTODOTODO

Wiki Link  
Video Link  
Youtube Link  

For more details, visit our team wiki. The system offers:
- Clear guidance through the MESA Receptor and ICD (Intracellular Domain) creation
- REST API, Python (PyPi) Package, Containerized Docker solutions and GitLab Repo-based deployment
- Great expandability through python
- A webapp and video guides
- Standardized Genbank file output
- SynBio specific and customizable sequence optimization

## Installation
This tool has been designed and developed with flexibility in mind. To achieve maximum compatibility, many different solutions exist.
All installation methods are listed below. For a quick deployment we recommend following the **pre-built docker images** installation instructions.

### Pre-Built Docker Image
#### Requirements
- docker: latest version. You can download and install this according to the official installation instructions: https://www.docker.com/products/docker-desktop/

#### Steps
1. You can preview pre-built docker images here[](https://hub.docker.com/repository/docker/aeneastews/mesa-designer/general) and find a preferred version (we recommend the latest).
Note: For most versions there are different containers available: **all**, **webapp** and **api**.
These can be denoted by their respective names ending in -all, -webapp and -api respectively. You can download (pull) the latest official version from dockerhub by running this command from a commandline or terminal while docker is running on your system:
```
docker pull aeneastews/mesa-designer:latest
```
2. This image can be started by running the following command from a commandline or terminal while docker is running on your system:
```
docker run -p 8501:8501 -p 8000:8000 aeneastews/mesa-designer:latest
```

### Custom Docker Image compilation
#### Requirements
- docker: latest version. You can download and install this according to the official installation instructions: https://www.docker.com/products/docker-desktop/

#### Steps
1. Download or clone the Dockerfile from the repository. This can be done by navigating to the Munich 2025 iGEM GitLab repo:
https://gitlab.igem.org/2025/software-tools/munich and downloading the [**Dockerfile**](https://gitlab.igem.org/2025/software-tools/munich/-/blob/main/Dockerfile) from the list of files.
2. Move the downloaded **Dockerfile** to a directory of your choice and ensure it is named "Dockerfile"
3. Start docker on your system if not already running. You can navigate to the same directory using a commandline or terminal and run the following command to compile the docker image:
```
docker build -t mesa-designer .
```
4. This image can be started by running the following command from a commandline or terminal while docker is running on your system:
```
docker run -p 8501:8501 -p 8000:8000 mesa-designer
```

### Local Installation
#### Requirements
- python: version=3.13.x. You can download and install this according to the official installation guides: https://www.python.org/downloads/. It may already be installed, you can check this by running `python --version` from a commandline or terminal.
- git: This is recommended for cloning the repository, though not strictly necessary. This will, however, enable easy updates to future versions.

#### Steps
1. Download or clone the repository:  
If you elected not to download git (not recommended), you have to download the zip archive from https://gitlab.igem.org/2025/software-tools/munich by clicking the Code dropdown.  
Assuming you downloaded git, run this command from a commandline or terminal:
```
git clone https://gitlab.igem.org/2025/software-tools/munich.git  
cd munich
```

2. Running the tool requires a few dependencies. To keep these separated, you should create a [python virtual environment](https://docs.python.org/3/library/venv.html). 
This can be done by running this command from a commandline or terminal. 
Note: you should do this inside the munich directory which was downloaded during the previous step:
```
python -m venv .venv
```
3. To install the required packages in the virtual environment, it has to be activated first.
```
On Windows run this command from the munich directory:
.\.venv\Scripts\activate

On Mac OS or Linux run this command from the munich directory:
source ./.venv/bin/activate
```
The required packages will be automatically installed and setup via this command on both systems:
```
pip install -r requirements.txt
```
4. The MESA Designer tool is based on validated experimental data from multiple databases.
This data can be automatically set up using the included setup.py script. Please run it using:
```
python setup.py
```
5. Since the system comes with two available services, you can start these individually to your liking.
Starting the webapp can be achieved by running the following command from a commandline or terminal:
```
On Windows run this command from the munich directory:
streamlit run .\app\main.py

On Mac OS or Linux run this command from the munich directory:
streamlit run ./app/main.py
```
6. Starting the API does not require the webapp to be running, however, if you wish to start both at the same time, we recommend doing so in a new commandline or terminal window.
From the munich directory run the following command:
```
uvicorn api.main:app
```

### Python (PyPi) Package
For easy integration with custom projects, the python package provides extended features, low-level control and complete documentation.
This method is intended to be used by developers seeking to integrate MESA functionality into their project.

#### Requirements
- python>=3.10
- biopython: compatible version with your own project

#### Steps
1. Having activated your desired installation python environment, run:
```
pip install mesa-designer
```


## Usage

### Webapp & API
Depending on your specific choice, you may have started different parts of the toolkit. This means that if you are running the local installation, you can choose which parts to run, and if you download a pre-built docker image, you will have different options available depending on your choice:
- -all: both
- -webapp: webapp
- -api: api  

Regardless of your choice, by default, the webapp can be accessed at: http://localhost:8501/ and the API at: http://localhost:8000/.
To get an overview of all API endpoints and their specific schemas, navigate to: http://localhost:8000/docs. This also provides you with the ability to try out all endpoints.

### Python Package
The python package can be used by importing it wherever necessary: 
```
import mesa_designer
```

Despite its extensive functionality, it is straightforward to use. Please check out the extensive annotation in the package's files and functions.
Some examples are provided here:

Getting an overview of preselected parts and components:
```python
from mesa_designer import ALL_DATA

print(ALL_DATA)
```

Building a new MESA chain, appending a few components and exporting it to an annotated genbank file:
```python
from mesa_designer.mesa import *

# create a new chain object
chain = MesaChain()

# add a custom binder design, this might be the result of a bindcraft https://github.com/martinpacesa/BindCraft run
chain.add_binder(sequence="BINDERSEQUENCE")

# add a default mesa tmd linker, this is a preconfigured linker, however the add_tmd_linker function and other linker functions provide extensive flexibility and customizability
chain.add_tmd_linker()

# the mutated TEVp N-Terminal end is a well-known and widely used N-Terminal end of a split TEV protease
chain.add_protease(protease_name="NTEVp_H75S")

# finally save the constructed mesa chain to an annotated genbank file
chain.save_genbank_file("mesa_chain.gb")
```

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started.
Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps
explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce
the likelihood that the changes inadvertently break something. Having instructions for running tests is especially
helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.  
streamlit   
biopython

Key Scientific References:

1. MESA Receptor Framework (Original):
   Daringer, N.M., Dudek, R.M., Schwarz, K.A., Leonard, J.N. (2014)
   "Modular Extracellular Sensor Architecture for Engineering Mammalian Cell-based Devices"
   ACS Synthetic Biology
   DOI: https://pubs.acs.org/doi/10.1021/sb400128g

2. MESA Applied to Human Cells:
   Schwarz, K.A., Daringer, N.M., Dolberg, T.B., Leonard, J.N. (2017)
   "Rewiring human cellular input-output using modular extracellular sensors"
   Nature Chemical Biology
   PubMed: https://pubmed.ncbi.nlm.nih.gov/27941759/

3. Elucidation and Refinement of Synthetic Receptor Mechanisms:
   Edelstein, H.I., Donahue, P.S., Muldoon, J.J., Kang, A.K., Dolberg, T.B., Battaglia, L.M., Allchin, E.R., Hong, M., Leonard, J.N. (2020)
   "Elucidation and refinement of synthetic receptor mechanisms"
   Synthetic Biology, 5(1), ysaa017
   DOI: https://doi.org/10.1093/synbio/ysaa017
   PubMed: https://pubmed.ncbi.nlm.nih.gov/33392392/

4. SAbDab - Structural Antibody Database:
   Dunbar, J., Krawczyk, K. et al. (2014)
   "SAbDab: the structural antibody database"
   Nucleic Acids Research, 42, D1140-D1146
   DOI: https://academic.oup.com/nar/article/42/D1/D1140/1044118

5. SKEMPI 2.0 - Protein Interaction Database:
   Jankauskaitė, J., Jiménez-García, B., Dapkūnas, J., Fernández-Recio, J., Moal, I.H. (2019)
   "SKEMPI 2.0: an updated benchmark of changes in protein–protein binding energy, kinetics and thermodynamics upon mutation"
   Bioinformatics, 35(3):462-469
   DOI: https://academic.oup.com/bioinformatics/article/35/3/462/5055583

Acknowledgments:
We thank all contributors and the scientific community for making this work possible.

## Design Outline
[Design Doc](https://docs.google.com/document/d/1ciPsgLo5JNp7wKqFREEnWSCiShK2ZrRBRCNdM3VBm_A/edit?tab=t.0)
