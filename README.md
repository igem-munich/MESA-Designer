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
now that the image is built, run it
```bash
docker run -p 8501:8501 mesa-designer 
```
### Pulling image from dockerhub
The MESA-Designer tool can now also be run as a docker container. Please find it [here](https://hub.docker.com/repository/docker/aeneastews/mesa-designer/general).

#### Running docker container
You can run the docker container using : `docker run -p 8501:8501 mesa-designer`

It can then be accessed from a webbrowser at localhost:8501

### Pre-Release Versions

##### v1.3-slim-beta
This version is the first to include an API. In order to access it you have to run the container using:   
`docker run -p 8501:8501 -p 8000:8000 mesa-designer`  
You can find the api documentation at http://localhost:8000/docs.


# old stuff



# Team Munich 2025 Software Tool

If your team competes in the [**Software & AI** village](https://villages.igem.org) or wants to
apply for the [**Best Software Tool** prize](https://competition.igem.org/judging/special-prizes), you **MUST** host all the
source code of your team's software tool in this repository, `main` branch. By the **Wiki Freeze**, a
[release](https://docs.gitlab.com/ee/user/project/releases/) will be automatically created as the judging artifact of
this software tool. You will be able to keep working on your software after the Grand Jamboree.

> If your team does not have any software tool, you can totally ignore this repository. If left unchanged, this
repository will be automatically deleted by the end of the season.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might
be unfamiliar with (for example your team wiki). A list of Features or a Background subsection can also be added here.
If there are alternatives to your project, this is a good place to list differentiating factors.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew.
However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing
specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a
specific context like a particular programming language version or operating system or has dependencies that have to be
installed manually, also add a Requirements subsection.

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
