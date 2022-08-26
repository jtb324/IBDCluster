---
layout: post
title:  "ChangeLog 8-25-22"
---

# Changelog - [1.1.0] - 2022-08-25:
---

## **Breaking:** Created a way for the user to create a custom plugins directory.

In previous iterations of the IBDCluster program, plugins were kept within the plugins directory found within the IBDCluster install directory. Users could add new plugins to this directory and then appropriately adjust the config.json file. The issue with this format is that users have to kep track of where this install directory is whenever they want to add new plugins. Now users can specify a custom directory anywhere they want and place the new plugins into this directory. The convenience of this approach is that the user can also use the stock plugins by just including those in the config.json. The user will have to change a few things to make the program work. These things are listed below. Since this is a breaking change the version is going to be updated to 1.2.0

**Creation of an IBDCLUSTER_MAIN_PLUGINS and IBDCLUSTER_CUSTOM_PLUGINS Environmental Variables:**
The user has to export two environmental variables in there .bashrc/.zshrc/(whatever shell you use config file). One variable has the absolute path to the stock plugins (IBDCLUSTER_MAIN_PLUGINS) and the other variable has the user custom defined directory (IBDCLUSTER_CUSTOM_PLUGINS). These export statements are shown below:

```bash
export IBDCLUSTER_CUSTOM_PLUGINS=...
export IBDCLUSTER_MAIN_PLUGINS=...
```

 IBDCluster will search for both of these variables in the environment namespace during runtime. A runtime error will not occur if neither directory is present. This happens because the method to get the environmental variable will just return None if it can't find it in the namespace. The error will occur when the program tries to match the directory with the config.json file. For the program to run the directories used in the config.json have to exist and the path has to be present in the shell configuration file.

 { : .optional}
 For just a stock install of IBCluster, the user just needs to install this line to the shell configuration file with the appropriate directories to the stock plugin directory:

 ```bash
export IBDCLUSTER_MAIN_PLUGINS={install_path}/{IBDCluster}/plugins/
 ```

**Adjusting the config.json file:**
An example config file is shown below. The first three plugins are stock plugins found in the plugin folder in the install directory. The fourth plugin is a custom plugin within the custom_plugins directory.

```json
{
    "plugins":["plugins.pvalues", "plugins.network_writer", "plugins.allpair_writer", "custom_plugins.phecode_pvalue_distributions"],
    "modules": [
        {
            "name":"pvalues"
        },
        {
            "name": "network_writer"
        },
        {
            "name": "allpair_writer"
        },
        {
            "name": "phecode_pvalue_distributions"
        }
    ]
}
```