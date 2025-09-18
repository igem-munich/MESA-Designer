# Team Munich 2025 Software Tool

## Description
The iGEM Munich 2025 MESA Designer Software Tool is an easy-to-use [MESA (Modular Extracellular Sensor Architecture for Engineering Mammalian Cell-based Devices)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4161666/) design tool.
The MESA Receptor Framework is a highly flexible and adaptable method for custom cellular circuit creation. The MESA Designer Tool provides a ... package for creation. Docker, Python package, webapp.
  
For more details, visit our team wiki. The system offers:
- Clear guidance through the MESA Receptor and ICD (Intracellular Domain) creation
- REST API, Python (PyPi) Package, Containerized Dockersolutions and GitLab Repo based deployment
- Great expanability through python
- A webapp and video guides
- Standardized Genbank file output
- SynBio specific and customizable sequence optimization

## Design Outline
[Design Doc](https://docs.google.com/document/d/1ciPsgLo5JNp7wKqFREEnWSCiShK2ZrRBRCNdM3VBm_A/edit?tab=t.0)

## Installation
This tool has been designed and development with flexibility in mind. To achieve maximum compatibility many different solutions exist.
All installation methods are listed below. For a quick deployment we recommend following the **pre-build docker images** installation instructions.

### Pre-Built Docker Image
#### Requirements
- docker: latest version. You can download and install this according to the official installation instructions: https://www.docker.com/products/docker-desktop/

#### Steps
1. You can preview pre-built docker images here [here](https://hub.docker.com/repository/docker/aeneastews/mesa-designer/general) and find a preferred version (we recommend the latest).
Note: For most versions there are different containers available: **all**, **webapp** and **api**.
These can be denoted by their respective names ending in -all, -webapp and -api respectively. You can download (pull) the latest official from dockerhub version by running this command from the a commandline or terminal while docker is running on your system:
```
docker pull aeneastews/mesa-designer:latest
```
2. This image can be started by running the following command from a commandline or terminal while docker is running on you system:
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
- git: This is recommended for cloning the repository though not strictly necessary. This will however enable easy updates to future versions.

#### Steps
1. Download or clone the repository:  
If you elected not to download git (not recommended), you have to download the zip archive from https://gitlab.igem.org/2025/software-tools/munich by clicking the Code dropdown.  
Assuming you downloaded git, run this command from a commandline or terminal:
```
git clone https://gitlab.igem.org/2025/software-tools/munich.git  
cd munich
```

2. Running the tool requires a few dependencies. To keep these separated you should create a [python virtual environment](https://docs.python.org/3/library/venv.html). 
This can be done by running this command from a commandline or terminal. 
Note: you should do this inside the munich directory which was downloaded during the previous step:
```
python -m venv .venv
```
3. To install the required packages in the virtual environment it has to be activated first.
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
This data can be automatically setup using the included setup.py script. Please run it using:
```
python setup.py
```
5. Since the system comes with two available services you can start these individually to your liking.
Starting the webapp can be achieved by running the following command from a commandline or terminal:
```
On Windows run this command from the munich directory:
streamlit run .\app\main.py

On Mac OS or Linux run this command from the munich directory:
streamlit run ./app/main.py
```

### Python (PyPi) Package
For easy integration with custom projects the python package provides extended features, low-level control and complete documention.
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
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of
usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably
include in the README.

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
