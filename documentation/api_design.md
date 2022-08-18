---
layout: page 
title: Api Design
parent: Customization
grand_parent: Documentation
nav_order: 1
---
# API Design:
---

The api is designed so that the user can access the following models as well as the output path for the program:

```python
@dataclass
class DataHolder:
    gene_name: str
    chromosome: int
    affected_inds: dict[str, list[str]]
    phenotype_prevalence: dict[str, float]
    phenotype_cols: list[str]
    ibd_program: str
    phenotype_description: None | dict[str, str] = None

@dataclass
class Pairs:
    """Object that has the data for each pair and contains a method
    to format the output string for the allpairs.txt file using this data"""

    pair_1: str
    phase_1: str
    pair_2: str
    phase_2: str
    chromosome: int
    segment_start: int
    segment_end: int
    length: float

    def form_id_str(self) -> str:
        """Method that will return the pair ids and phases in a formated string"""
        return f"{self.pair_1}\t{self.pair_2}\t{self.phase_1}\t{self.phase_2}"

    def form_segment_info_str(self) -> str:
        """Method that will return a string with all of the segment info"""
        return f"{self.segment_start}\t{self.segment_end}\t{self.length}\n"

@dataclass
class Network:
    """This class is going to be responsible for the clustering of each network"""

    gene_name: str
    gene_chr: str
    network_id: int
    pairs: list[Pairs] = field(default_factory=list)
    iids: set[str] = field(default_factory=set)
    haplotypes: set[str] = field(default_factory=set)
    connections: dict[str, int] = field(default_factory=dict)
    pvalues: None | str = None
```

The api will supply three values: "data", "network", "output" which are represented by the DataHolder model, the Network model, and the gene output path respectively. These attributes and methods of each model are described below.

## DataHolder: 
---
The attributes for the DataHolder model are described below:

* <span style="color: #F0FF00">**gene_name**</span> - string that contains the gene name that was passed to the program in the gene_info file

* <span style="color: #F0FF00">**chromosome**</span> - integer number for the chromosome of interest

* <span style="color: #F0FF00">**affected_inds**</span> - dictionary where the keys are a string for each PheCode and the values are a list of IDs that are affected by the phenotype.

* <span style="color: #F0FF00">**phenotype_prevalence**</span> - dictionary where the keys are a PheCode string and the values are the percentages of individuals that are affected by the phenotype in the population you provided.

* <span style="color: #F0FF00">**phenotype_cols**</span> - list of PheCode strings in the provided PheCode matrix.

* <span style="color: #F0FF00">**ibd_program**</span> - this is a string that gives the name of the ibd program being used.

* <span style="color: #F0FF00">**phenotype_description**</span> - This attribute is either none if the user didn't provide a file that has descriptions of each PheCode or it is a dictionary where the keys are PheCode strings and the values are the description of each PheCode.

## Pairs:
---
The attributes for the Pairs model are described below:

* <span style="color: #F0FF00">**pair_1**</span> - string that contains the ID of the first person in the pair that shares an IBD segment.

* <span style="color: #F0FF00">**phase_1**</span> - string that contains ID and the phase number of pair_1's haplotype which contains the shared IBD segment. This values is the pair ID concatenated with the phase value. Ex: ID1 becomes ID1.1 for the first haplotype (Different ibd detection programs use different numbers to represent haplotypes. The example above is based on hap-IBD).

* <span style="color: #F0FF00">**pair_2**</span> - string that contains the ID of the second person in the pair that shares an IBD segment.

* <span style="color: #F0FF00">**phase_2**</span> - string that contains the ID and the phase number of pair_2's haplotype which contains the shared IBD segment. This values is the pair ID concatenated with the phase value. Ex: ID2 becomes ID2.2 for the second haplotype (Different ibd detection programs use different numbers to represent haplotypes. The example above is based on hap-IBD).

* <span style="color: #F0FF00">**chromosome**</span> - integer number of the chromosome that the shared segment is on.

* <span style="color: #F0FF00">**segment_start**</span> - integer that has the start base position of the shared segment.

* <span style="color: #F0FF00">**segment_end**</span> - integer that has the end base position of the shared segment.

* <span style="color: #F0FF00">**length**</span> - float value that has the total length of the segments in centimorgans.

The methods for the Pairs model are described below:
* <span style="color: #00FFE0">**form_id_str**</span> - method that creates a string with tab separated values that has the pair_1 ID string, the pair_2 ID string, the phase_1 string, and the phase_2 string, respectively.

* <span style="color: #00FFE0">**form_segment_info_str**</span> - method that creates a string with tab separated values that has the segment start position, the segment end position, and the length of the segment, respectively.

## Network:
---
The attributes for the Network model are described below:

* <span style="color: #F0FF00">**gene_name**</span> - string that contains the gene name that was passed to the program in the gene_info file.

* <span style="color: #F0FF00">**gene_chr**</span> - integer that has the chromosome number that the gene is on.

* <span style="color: #F0FF00">**network_id**</span> - integer that has the id for the network in the analysis. 

* <span style="color: #F0FF00">**pairs**</span> - list of Pair objects that have information about the ibd segment that each id pair shares.

* <span style="color: #F0FF00">**iids**</span> - a set that has all the unique iid strings that are in the network.

* <span style="color: #F0FF00">**haplotypes**</span> - a set that has all the unique haplotypes in the network. The unique haplotypes are just the ID with the phase information appended to it. Ex: ID1 becomes ID1.1. This attribute allows us to differentiate between which haplotype is in the network. 

* <span style="color: #F0FF00">**connections**</span> - a dictionary where the keys are the individual grid strings and the values are the number of other individuals that the ID is connected to.

* <span style="color: #F0FF00">**pvalues**</span> - This value is either None or it is a string of tab separated values where the first three values are the most significant p-value, the corresponding PheCode, and then the corresponding PheCode description respectively. The rest of the string are tab separated values for the p-value of each PheCode and how many individuals in the network are affected by the Phenotype. The default value for this attribute is None and the pvalue.py plugin actually adds the string to this attribute. The user could substitute there own statistically plugin if they wished and use this attribute, it just needs to meet the plugin requirements which are defined in the "Custom Pluging Design" section.

## Output Path:
---
The api also supplies the gene output path. This output path is formed by affixing the name of the gene of interest with the output path provided by the user to the CLI. The plugin needs to create the path to write information to.

for example if you were running the gene rbm20, then this value would be:

```
/output_path_provided_by_user/rbm20/
```