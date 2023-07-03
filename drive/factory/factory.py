from typing import Any, Callable, Protocol


class PluginNotFound(Exception):
    """
    Error that is raised if the user tries to load a plugin
    that is not there
    """

    def __init__(self, plugin_type: str) -> None:
        super().__init__(
            f"Plugin, {plugin_type} not found in the plugins folder. Make sure you loaded the plugin into the json file."  # noqa: E501
        )


class AnalysisObj(Protocol):
    """Interface defining of an analysis object"""

    def analyze(self, **kwargs) -> Any:
        """
        Method that will analyze the inputs according to the purpose
        of the plugin
        """


analyze_obj_creation_funcs: dict[str, Callable[..., AnalysisObj]] = {}


def register(plugin_name: str, creation_func: Callable[..., AnalysisObj]) -> None:
    """registers the AnalysisObj plugin"""
    analyze_obj_creation_funcs[plugin_name] = creation_func


def unregister(plugin_name: str) -> None:
    """function that will unregister the plugin

    Parameters
    ----------
    plugin_name : str
        name of the plugin to be loaded. This name will be in json file
    """
    analyze_obj_creation_funcs.pop(plugin_name, None)


def create(arguments: dict[str, Any]) -> AnalysisObj:
    args_copy = arguments.copy()

    plugin_type: str = args_copy.pop("name")

    try:
        creation_func = analyze_obj_creation_funcs[plugin_type]
        return creation_func(**args_copy)
    except KeyError:
        raise PluginNotFound(plugin_type)  # pylint: disable=raise-missing-from
