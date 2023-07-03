"""File used to load in the different plugins"""

import importlib
from typing import Protocol


class PluginInterface(Protocol):
    """Interface that will define how a plugin looks like"""

    @staticmethod
    def initialize() -> None:
        """Method that will initialize the plugin"""


def import_module(name: str) -> PluginInterface:
    """Import the plugin so it can be initialized

    Parameters
    ----------
    name : str
        name of the plugin

    Returns
    -------
    PluginInterface
    """
    return importlib.import_module(name)  # type: ignore


def load_plugins(plugins: list[str]) -> None:
    """Calls the initialize method for each plugin"""
    for plugin_name in plugins:
        plugin = import_module(plugin_name)

        plugin.initialize()
