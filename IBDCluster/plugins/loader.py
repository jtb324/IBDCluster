"""File used to load in the different plugins"""

import importlib
from typing import Protocol, List


class PluginInterface(Protocol):
    """Interface that will define how a plugin looks like"""

    @staticmethod
    def initialize() -> None:
        """Method that will initial the plugin"""


def import_module(name: str) -> PluginInterface:
    return importlib.import_module(name)  # type: ignore


def load_plugins(plugins: List[str]) -> None:
    """Calls the initialize method for each plugin"""
    for plugin_name in plugins:
        plugin = import_module(plugin_name)

        plugin.initialize()
