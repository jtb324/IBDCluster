Custom Plugins
==============
Users are free to design their own plugins to fit their own needs, but it has to conform to the following design.

Plugin API Structure
--------------------

DRIVE uses an object oriented approach for the plugin architecture. There are ways to achieve this architecture without classes but the developer choose this structure so if you are not familiar with this programming paradigm in python you might want to read up on it before continuing. 


Required methods/functions/attributes:
--------------------------------------

DRIVE requires the plugin class to have one method called "analyze". This method is called by the main DRIVE code. This method should only accept the arguments "self" and "\*\*kwargs and should not return a value. This kwargs argument will contain the Data container that has the networks identified from the clustering analysis. This container is described in further detail in a different section. This is the only mandatory method. Users can add other methods but they need to be called within this analyze method because DRIVE will only call the analyze method. This is the only mandatory method. Users can add other methods but they need to be called within this analyze method because DRIVE will only call the analyze method. 

.. note:: 

    Although the analyze function cannot return a value. It can modify the data container and this modification would be used in any plugins that follow this one.


DRIVE also requires the plugin class to have an attribute called "name". This attribute is used for logging messages so DRIVE will crash if it is not present.

This plugin file is required to have one additional function called "initialize" which accepts no arguments and returns nothing. This function is used by DRIVE to run the code dynamically. 

A template for the plugin class is shown below. The example uses dataclasses but a traditional class will work as well

..  code-block:: python

    from dataclasses import dataclass

    # The way that DRIVE implements the plugin architecture requires that an object 
    # oriented approach be used. We are showing data classes in this example but 
    # tradition python classes will work as well

    @dataclass
    class PluginName:

        # DRIVE requires the plugin to have an attribute called name. This attribute is 
        # used for logging
        name: str = "Name of the plugin"

        def analyze(self, **kwargs) -> None:
            """This is the main function that DRIVE will use to interact with 
            the plugin. This function is required and has to be called analyze"""

            ...


    # The pluging file needs to have a function at the end called initialize. This 
    # function is used by the plugin factory to dynamically run the code at runtime.
    def initialize() -> None:
        factory_register("plugin_name from the file name", PluginName)

Where the plugins have to be place:
-----------------------------------
DRIVE expects the code for the plugins to be located in a module within the DRIVE source code directory. Within this source code there is a subdirectory call plugins. This directory is where DRIVE will look for the plugin code.

.. attention::
    At the moment this is the only directory where DRIVE looks for plugins. In the future, there may be added support for the user to specify another directory where plugins may be located.