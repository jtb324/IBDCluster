---
layout: page 
title: Stock Plugins
parent: Customization
grand_parent: Documentation
nav_order: 2
---
# Stock plugins:
---
IBDCluster comes with three stock plugins. Each of these plugins access the api described by the "api_design" section of the documentation. These plugins are responsible for the calculation of the different pvalues for each phenotype, writing out the information for each network to a file, and writing out the information for each of the pairs to the file. All of these plugins are found within the plugin directory within the ibd source code directory. This directory can be found within where you installed the program. Once inside the install directory if you just run the following two commands then you should see the plugins subdirectory:

```bash
cd IBDCluster/
```

```bash
ls
```
This directory is designed as a python module so it has an \_\_init\_\_.py file so that the python interpreter can import the code. The other three files in these directory are the plugins.

## Overview of the Plugin design:
---

### Basic structure of the plugin class:

The stock plugins use an object oriented design meaning each plugin is a class that has a method that gets called by the IBDCluster program. An example of the basic required plugin file structure is shown below and we'll discuss each part.

```python
from dataclasses import dataclass
from factory import factory_register

@dataclass
class Plugin_Class:
    """Example of the base requirements for a plugin class"""

    name: str = ""

    def analyze(self, **kwargs) -> None:
        """the analysis method is what gets called by the IBDCluster program"""

def initialize() -> None:
    factory_register("pvalues", Pvalues)
```

***Imports:***<br>
The example plugin here imports dataclass and a method called factory_register. If you are unfamiliar with dataclasses you can read about it here, ([Documentation link here](https://docs.python.org/3/library/dataclasses.html)), but they are basically just newer abstractions of the standard python class (with a bit of subjective pros/cons). You can use a normal class structure, you just need to create a \_\_init\_\_ function with the name attribute. The plugin also imports the factory_register function. We're not going to describe this plugin in detail because it just abstracts away some of the core functionality of the IBDCluster program. This function comes from a program module called factory (This module is not a standard python library so please don't try to find it). This function is required to have your module properly loaded into the program. 

***Class Structure:***<br>
The plugin class only has one required attribute and this the plugin's name. This name is used in a logging statement by the program (Therefore the program will give you an error if this is not there). The plugin also has to have one method named "analyze" which only takes key word arguments. This method also doesn't return anything. The three keyword arguments that the program passes to this method, will be "data", "network", "output". These arguments are what the plugins access. 

***Initialize function:***<br>
Each plugin file also has a function called initialize. This function doesn't take any inputs and returns nothing. The purpose of this function is to called the factory_register function. It passes the name of the plugin as the first argument and the plugin Class as the second argument.

{: .warning }
***Warning:***<br>
The name argument that is passed to the factory_register has to correspond to the name in the methods section of the config.json file. If it is not the same then the program will give an error when trying to load in the plugin. 

{: .optional }
***Note on Plugin order:*** <br>
The order of these plugins matters (To some degree) and that will be described below.

## Overview of each Plugin:
---

### pvalues.py:
This plugin is responsible for determining the enrichment p-values from the binomial test for each phenotype for each network. To determine the p-values, the program relies on the binomtest function from the scipy.stats library. The plugin also determines the minimum p-value found for each network as well as the phenotype that the p-value represents. The plugins creates a string with all of this information. A example string of what will be formed for two phenotypes is shown below:

```python
f"{minimum pvalue}/t{correspond phenotype}\t{phenotype description}\t{number of people affected by phenotype A}\t{phenotype A p-value}\t{number of people affected by phenotype B}\t{phenotype B p-value}"
```

This string will be saved in the network.pvalues attribute of the network models.

### network_writer.py:
Function that writes all of the information to the networks.txt output file. The description of this output file can be found in the "outputs" section of the documentation. 

### allpair_writer.py:
Function that writes all of the information to the allpairs.txt output file. The description of this output file can be found in the "outputs" section of the documentation. 