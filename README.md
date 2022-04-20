[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# IBDCluster v0.9.1:

## Purpose of the project: 
___
This project is a cli tool that clusters shared ibd segments within biobanks around a gene of interest. These network are then analyzed to determine how many individuals within a network are affected by a phenotype of interest.

## General PipeLine:
___
```mermaid
flowchart LR
    A(IBD Information) --> B(Identified Networks) --> C(Binomial test for enrichment of Phenotypes)
```

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
* You can find all the optional parameters by running:
```
IBDCluster.py --help
```
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
│   ├── tests
│   │   ├── test_data
│   │   ├── test_integration

```
## Comments about models:
___
* Classes for the cluster_class.py:

```mermaid
classDiagram
    class Cluster {
        ibd_file: str
        ibd_program: str
        indices: models.FileInfo
        count: int=0
        ibd_df: pd.DataFrame=pd.DataFrame
        network_id: str=1
        inds_in_network: Set[str]=set
        network_list: List[Network]=list
    }
    class Network {
        gene_name: str
        gene_chr: str
        network_id: int
        pairs: List[Pairs]=list
        iids: Set[str]=set
        haplotypes: Set[str]=set
        +filter_for_seed(ibd_df: pd.DataFrame, ind_seed: List[str], indices: FileInfo, exclusion: Set[str]=None) -> pd.DataFrame
        #determine_pairs(ibd_row: pd.Series, indices: FileInfo) -> Pairs
        +gather_grids(dataframe: pd.DataFrame, pair_1_indx: int, pair_2_indx: int) -> Set[str]
        +update(ibd_df: pd.DataFrame, indices: FileInfo) -> None
    }
    class FileInfo {
        <<interface>>
        id1_indx: int
        ind1_with_phase: int
        id2_indx: int
        ind2_with_phase: int
        chr_indx: int
        str_indx: int
        end_indx: int
        +set_program_indices(program_name: str) -> None
    }
    Cluster o-- Network
```

## Entity relationships:
___
```mermaid
erDiagram
    NETWORK }|--|{ PAIRS : contains
    NETWORK {
        string gene_name
        string chromosome
        int network_id
    }
    NETWORK }|--|{ IIDS : contains
    NETWORK }|--|{ HAPLOTYPES : contains
    PAIRS {
       string pair_1_id
       string pair_1_phase 
       string pair_2_id
       string pair_2_phase 
       int chromosome_number
       int segment_start 
       int segment_end
       float length 
       series affected_statuses 
    }
    IIDS {
        string Individual-ids
    }
    HAPLOTYPES {
        string haplotype-phase
    }
```
