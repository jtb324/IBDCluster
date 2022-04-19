from typing import Protocol, Dict, Callable, Any


class PluginNotFound(Exception):
    """
    Error that is raised if the user tries to load a plugin
    that is not there
    """

    def __init__(self, plugin_type: str) -> None:
        super().__init__(
            f"Plugin, {plugin_type} not found in the plugins folder. Make sure you loaded the plugin into the json file."
        )


class AnalysisObj(Protocol):
    """Basic representation of an Analysis object"""

    def analyze(self) -> None:
        """
        Method that will analyze the inputs according to the purpose
        of the plugin
        """


analyze_obj_creation_funcs: Dict[str, Callable[..., AnalysisObj]] = {}


def register(plugin_name: str, creation_func: Callable[..., AnalysisObj]) -> None:
    """registers the AnalysisObj plugin"""
    analyze_obj_creation_funcs[plugin_name] = creation_func


def unregister(plugin_name: str) -> None:
    """function that will unregister the plugin"""
    analyze_obj_creation_funcs.pop(plugin_name, None)


def create(arguments: Dict[str, Any]) -> AnalysisObj:

    args_copy = arguments.copy()

    plugin_type: str = args_copy.pop("name")

    try:
        creation_func = analyze_obj_creation_funcs[plugin_type]
        return creation_func(**args_copy)
    except KeyError:
        raise PluginNotFound(plugin_type)
