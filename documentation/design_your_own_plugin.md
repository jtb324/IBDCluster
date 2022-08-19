---
layout: page 
title: Custom Pluging Design
parent: Customization
grand_parent: Documentation
nav_order: 3
---
# Designing your own plugin:
---

When you install the program you need to add two environmental variables to your .bashrc file (or .zshrc or whatever shell you use). These two variables are the IBDCLUSTER_MAIN_PLUGINS and the IBDCLUSTER_CUSTOM_PLUGINS. The IBDCLUSTER_MAIN_PLUGINS has filepath to the directory with the stock plugins and then the IBDCLUSTER_CUSTOM_PLUGINS has the filepath to the directory with the users custom plugins. These two directories can not have the same name or the program will not work. 