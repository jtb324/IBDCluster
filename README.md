# IBDCluster v0.9.1:

## Purpose of the project: 
___
This project is a cli tool that clusters shared ibd segments within biobanks around a gene of interest. These network are then analyzed to determine how many individuals within a network are affected by a phenotype of interest.

## installing:
___
***Cloning from github and modify permissions:***
1. Clone the project into the appropriate directory using git clone.
2. cd into the IBDCluster directory
```
cd IBDCluster
```
2. run the following command to set the right permissions on the IBDCluster.py file
```
chmod +x IBDCluster/IBDCluster.py
```
***Installing dependencies:***
Next install all the necessary dependencies. The steps for this vary depending on what package manager you are using.

*If using conda:*
1. There is a environment.yml file in the main IBDCluster directory. Run the following command and it will create an environment called IBDCluster

```
conda env create -f environment.yml
```

2. You can now activate the environment by calling:

```
conda activate IBDCluster
```

*If using mamba:*
1. This is the same as the conda section except use the command
```
mamba env create -f environment.yml
```
2. You can activate this environment using:
```
conda activate IBDCluster
```

*If using Poetry*
1. The requirements for a poetry project are also in the IBDCluster directory. Ideally you need to activate some type of virtual environment first. This environment can be either a conda environment or a virtualenv. Once this environment is activated you can call:

```
poetry install
```

2. At this point all necessary dependencies should be installed.

* if you wish to find more information about the project you can find the documentation here: https://python-poetry.org/

***Adding IBDCluster to the users $PATH:***
To be able to run the IBDCluster program without having to be in the source code directory, you should add the IBDCluster.py file to your path.

1. In your .bashrc file or .zshrc add the line :
```
export PATH="{Path to the directory that the program was cloned into}/IBDCluster/IBDCluster:$PATH"
```
2. run this line:
```
source .bashrc
```
or
```
source .zshrc
```
This will allow you to run the code by just typing IBDCluster.py from any directory.

***Running IBDCluster***

## Running the code:
___
*

## Technical Details of the project:
___
* This part is mainly for keeping track of the directory structure.

## Project Structure:
___
```
├── IBDCluster
│   ├── analysis
│   │   ├── main.py
│   │   ├── percentages.py
│   ├── callbacks
│   │   ├── check_inputs.py
│   ├── models
│   │   ├── cluster_class.py
│   │   ├── indices.py
│   │   ├── pairs.py
│   │   ├── writers.py
│   ├── log
│   │   ├── logger.py
│   ├── cluster
│   │   ├── main.py
│   ├── IBDCluster.py
├── .env
├── environment.yml
├── .gitignore
├── poetry.lock
├── pyproject.toml
├── README.md
├── requirements.txt
```