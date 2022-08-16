---
layout: page
title: Github Installation
parent: Installation
grand_parent: Documentation
---
# Using Github to install IBDCluster:
---

This installation method assumes that you are familiar with Git and Github, the commandline, and python's Anaconda package manager and that these programs are install/can be install on whatever computing environment you are using. If you wish to read the documentation for each of these then they will be listed below:

* **Git:** [Git Website](https://git-scm.com/)

* **Github:** [Github Website](https://github.com/)

* **Commandline Interface:** [This is probably overkill but here is a very indepth CLI tutorial](https://www.learnenough.com/command-line-tutorial)

* **Anaconda** [Anaconda Website](https://www.anaconda.com/)

*Optional Installation Dependency:*
You can also use poetry to install the program. Poetry is a python package manager (another alternative to Pip and Conda and all the other package manages) that has good dependency resolution to create a reproducible environment. You can read more about the project here [Poetry documentation](https://python-poetry.org/) and the steps to install it are described here [Poetry Installation](https://python-poetry.org/docs/#installation). Poetry is the current recommended way to install the program but it relies on you having the necessary dependencies to install the Poetry program into whatever system you are using. People trying to run this on a personal machine (probably unlikely) should have no difficultly installing poetry but those running IBDCluster on a server (Which considering its made for BioBank data will be most people) may run into permission errors. If this is the case then just skip to the section titled "Creating a Conda Virtual Environment" and start the installation there.

