---
layout: page 
title: Customization
parent: Documentation
has_children: true
nav_order: 3
---
# Customizing the IBDCluster network analysis:
---
This program was built so that individuals can bring in custom analysis scripts using plugins. This allows flexibility for what enrichment test the user wishes to use as well as what information they would like to keep to output files. Three defaults plugins are provided and are used: pvalues.py, network_writer.py, allpair_writer.py. These plugins will be discussed in more detail later. The following sections will discuss the api design, the stock plugins, and how the user can design their own plugin.
