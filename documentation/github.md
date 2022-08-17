---
layout: page
title: Github Installation
parent: Installation
grand_parent: Documentation
nav_order: 1
---
# Using Github to install IBDCluster:
---

This installation method assumes that you are familiar with Git and Github, the commandline, and python's Anaconda package manager and that these programs are install/can be install on whatever computing environment you are using. You will have to use all of these tools so you will need to be familiar enough with each one to run the example commands. If you wish to read the documentation for each of these then they will be listed below:

* **Git:** [Git Website](https://git-scm.com/)

* **Github:** [Github Website](https://github.com/)

* **Commandline Interface:** [This is probably overkill but here is a very indepth CLI tutorial](https://www.learnenough.com/command-line-tutorial)

* **Anaconda** [Anaconda Website](https://www.anaconda.com/)

{: .optional }
***Optional Installation Dependency:*** <br> 
You can also use poetry to install the program. Poetry is a python package manager (another alternative to Pip and Conda and all the other package manages) that has good dependency resolution to create a reproducible environment. You can read more about the project here [Poetry documentation](https://python-poetry.org/) and the steps to install it are described here [Poetry Installation](https://python-poetry.org/docs/#installation). Poetry is the current recommended way to install the program but it relies on you having the necessary dependencies to install the Poetry program into whatever system you are using. People trying to run this on a personal machine (probably unlikely) should have no difficultly installing poetry but those running IBDCluster on a server (Which considering its made for BioBank data will be most people) may run into permission errors. If you are running into permission errors than the documentation will indicate where your installation instructions are different.

## Steps to installing IBDCluster:
The installation process can be broken into 4 steps. These are listed below and will be explained in further detail:

1. clone the Github repository to your local environment.
2. create and activate a conda virtual environment.
3. install all the dependencies using Conda or Poetry.
4. Add the program to you PATH so that you can call the program

## Step 1: Clone the Github repository:
You can clone the Github repository into your local environment using the command shown below:

```bash
git clone https://github.com/jtb324/IBDCluster.git
```

You should now have a directory called IBDCluster. You can check if this exists using the command:

```
ls IBDCluster/
```
If you see a directory file tree then the program cloned correctly. If you receive an error saying that the directory does not exist, then you will have to debug the error to move onto step 2.

## Step 2: create and activate a conda virtual environment:
To create a conda virtual environment you just need to run the command below, replacing the name with whatever name you want to give the environment.

```bash
conda create -n <environment-name> python=3.10
```

{: .warning }
IBDCluster requires the python version to be >= 3.10 because it uses features that were introduced in that release. If you try to run it in python=3.9 it will give many errors.

{: .alternative }
If you are not using poetry to install dependencies then you can skip step 2 and just move onto step 3.

## Step 3: install all the dependencies using Conda or Poetry:
---


## Step 4: Add the program to you PATH so that you can call the program:
---
