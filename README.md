# Team Munich 2025 Software Tool

## Description
The iGEM Munich 2025 MESA Designer Software Tool is an easy-to-use [MESA (Modular Extracellular Sensor Architecture for Engineering Mammalian Cell-based Devices)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4161666/) design tool.
The MESA Receptor Framework is a highly flexible and adaptable method for custom cellular circuit creation. The MESA Designer Tool provides a package for easy and fast MESA Receptor creation. Read more on our iGEM Wiki.

Check out our [Wiki](https://2025.igem.wiki/munich/)

For more details, visit our team wiki. The system offers:
- Clear guidance through the MESA Receptor and ICD (Intracellular Domain) creation
- REST API, Python (PyPi) Package, Containerized Docker solutions and GitLab Repo-based deployment
- Great expandability through python
- A webapp and video guides
- Standardized Genbank file output
- SynBio specific and customizable sequence optimization

## Installation Guides

Various installation guides exist, each with upsides and downsides. We will go through common deployment options in the following sections. If you wish to get started using MESA Designer as quickly as possible and do not care about modifying the source code or programmatic integration, we recommend using a **Pre-Built Docker Image**.

### Containerized Deployment

Docker deployment ensures consistency across environments and greatly simplifies installation. This also enables easy integration with existing environments via the use of [docker compose](https://docs.docker.com/compose/intro/features-uses/) scripts. This modular, containerized approach minimizes resource usage - a research group might deploy only the API on their compute cluster while running the web interface on a local workstation.

**Pre-Built Docker Images (Recommended for Beginners):**

#### Requirements
- Docker: latest version. You can download and install this according to the [official installation instructions](https://www.docker.com/products/docker-desktop/).

#### Steps
1. You can view pre-built docker images [on dockerhub](https://hub.docker.com/repository/docker/aeneastews/mesa-designer/general) and find a preferred version (we recommend the latest).  
Note: For most versions there are different containers available: **all**, **webapp** and **api**. These can be identified by their respective names ending in -all, -webapp and -api respectively. You can download (pull) the latest official version from dockerhub by running this command from a commandline or terminal while docker is running on your system:
```bash
docker pull aeneastews/mesa-designer:latest
```

2. The image can be started by running the following command from a commandline or terminal while docker is running on your system:
```bash
docker run -p 8501:8501 -p 8000:8000 aeneastews/mesa-designer:latest
```

3. You can now access the web application at [http://localhost:8501/](http://localhost:8501) and the API at [http://localhost:8000/](http://localhost:8000/)


**Production Deployment / Custom Build:** 

For institutional deployment, we provide easily adaptable dockerfiles for orchestrated multi-container deployment with reverse proxies, SSL termination, and automated backups. The software runs efficiently on standard cloud platforms such as [AWS EC2](https://aws.amazon.com/ec2/), [Google Cloud Platform](https://cloud.google.com/) or [Hetzner Cloud](https://www.hetzner.com/). The environments are also compatible with standard academic computing environments including everything from simple machines to [SLURM](https://en.wikipedia.org/wiki/Slurm_Workload_Manager) clusters.

#### Requirements
- Docker: latest version. You can download and install this according to the [official installation instructions](https://www.docker.com/products/docker-desktop/).

#### Steps
1. Download or clone the Dockerfile from the repository. This can be done by navigating to the [Munich 2025 iGEM GitLab repo](https://gitlab.igem.org/2025/software-tools/munich) and downloading the [**Dockerfile**](https://gitlab.igem.org/2025/software-tools/munich/-/blob/main/Dockerfile) from the list of files. You can then adapt the file to your specific needs or use the default settings for an easy deployment. However, in the latter case we recommend using a pre-built image.

2. Move the downloaded **Dockerfile** to a directory of your choice and ensure it is named **"Dockerfile"**

3. Start docker on your system if not already running. You can navigate to the same directory using a commandline or terminal and run the following command to compile the docker image:
```bash
docker build -t mesa-designer .
```

4. This image can be started by running the following command from a commandline or terminal while docker is running on your system:
```bash
docker run -p 8501:8501 -p 8000:8000 mesa-designer
```

5. You can now access the web application at [http://localhost:8501/](http://localhost:8501) and the API at [http://localhost:8000/](http://localhost:8000/)


### Local Installation

For researchers requiring local installation due to data sensitivity or customization needs, we provide the option of local deployment. This approach enables complete access, customizability and data privacy by removing the need for cross-device communication thus enabling increased security and protecting user's data.

#### Requirements
- Python: version=3.13.x. You can download and install this according to the [official installation instructions](https://www.python.org/downloads/). It may already be installed, you can check this by running `python --version` or `python3 --version` from a commandline or terminal. If the result is anything other than `Python 3.13.x # Note: x represents a flexible number` then you likely do not have the correct version of Python installed.
- Git: This is recommended for cloning the repository, though not strictly necessary. This will, however, enable easy updates to future versions. You can download it according to the [official installation instructions](https://git-scm.com/downloads). This may already be installed on your system. You can check this by running `git --version` from a commandline or terminal. If the result is anything other than `git version x.x.x` you will probably have to install it.

#### Steps
1. Download or clone the repository:  
If you elected not to download git (not recommended), you have to download the zip archive from [Munich 2025 iGEM GitLab repo](https://gitlab.igem.org/2025/software-tools/munich) by clicking the Code dropdown and selecting the zip-archive option.  
Assuming you downloaded git (recommended), run these commands from a commandline or terminal:
```bash
git clone https://gitlab.igem.org/2025/software-tools/munich.git  
cd munich
```

2. Running the tool requires a few dependencies. To keep these separated, you should create a [python virtual environment](https://docs.python.org/3/library/venv.html). 
This can be done by running this command from a commandline or terminal.  
Note: You should do this inside the munich directory which was downloaded during the previous step:
```bash
python -m venv .venv
```

3. To install the required packages in the virtual environment, it has to be activated first.  
```bash
# On Windows run this command from the munich directory:
.\.venv\Scripts\activate

# On Mac OS or Linux run this command from the munich directory:
source ./.venv/bin/activate
```
The required packages will be automatically installed and setup via this command on both systems:
```bash
pip install -r requirements.txt
```

4. The MESA Designer tool is based on validated experimental data from multiple databases. This data can be automatically set up using the included setup.py script. Please run it using:
```bash
python setup.py
```

5. Since the system comes with two available services, you can start these individually to your liking. Starting the webapp can be achieved by running the following command from a commandline or terminal:
```bash
# On Windows run this command from the munich directory:
streamlit run .\app\main.py

# On Mac OS or Linux run this command from the munich directory:
streamlit run ./app/main.py
```

6. Starting the API does not require the webapp to be running. However, if you wish to start both at the same time, we recommend doing so in a new commandline or terminal window.
From the munich directory run the following command in a new commandline or terminal window:
```bash
uvicorn api.main:app
``` 

7. You can now access the web application at [http://localhost:8501/](http://localhost:8501) and the API at [http://localhost:8000/](http://localhost:8000/)

### Python Package Index (PyPi) Package
For easy integration with custom projects, the python package provides extended features, low-level control and complete documentation. This method is intended to be used by developers seeking to integrate MESA functionality into their project or extend the MESA Designer Toolkit.

#### Requirements
- Python: version>=3.10. You can download and install this according to the [official installation instructions](https://www.python.org/downloads/). It may already be installed, you can check this by running `python --version` or `python3 --version` from a commandline or terminal. If the result is anything other than `Python 3.13.x # Note: x represents a flexible number` then you likely do not have the correct version of Python installed.
- biopython: compatible version with your own project. You can check it out on and install it from [PyPi](https://pypi.org/project/biopython/).

#### Steps
1. Having activated your desired installation python environment, run:
```bash
pip install mesa-designer
```

## Usage

For usage instructions, please checkout the "MESA Designer Usage Guide" in the guides directory and our videos on the igem wiki.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started.
Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps
explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce
the likelihood that the changes inadvertently break something. Having instructions for running tests is especially
helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
We want to thank the creators of [streamlit](https://streamlit.io/) for making the front-end development of modern applications straight forward, while allowing immense customizability.
   
We want to thank the creators of [biopython](https://biopython.org/) for enabling everyone to interact with many biological formats with ease, allowing developers to focus on functionality.

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

We thank all contributors and the scientific community for making this work possible.

## Design Outline
[Design Doc](https://docs.google.com/document/d/1ciPsgLo5JNp7wKqFREEnWSCiShK2ZrRBRCNdM3VBm_A/edit?tab=t.0)
