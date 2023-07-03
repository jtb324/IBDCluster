from dataclasses import dataclass, field
from typing import Any, Dict, List, Protocol, Set, TypeVar, Union

T = TypeVar("T", bound="Network")


class Network_Interface(Protocol):
    clst_id: float  # I don't like this attribute being a float but for now it has to remain this way for backwards compatibility
    true_positive_count: int
    true_positive_percent: float
    false_negative_edges: List[int]
    false_negative_count: int
    members: Set[int]
    haplotypes: List[int]
    min_pvalue_str: str = ""
    pvalues: Dict[str, str] = field(default_factory=dict)

    def print_members_list(self) -> str:
        """Returns a string that has all of the members ids separated by space

        Returns
        -------
        str
            returns a string where the members list attribute
            is formatted as a string for the output file. Individuals strings are joined by comma.
        """
        ...


@dataclass
class Network:
    clst_id: float
    true_positive_count: int
    true_positive_percent: float
    false_negative_edges: List[int]
    false_negative_count: int
    members: Union[Set[int], Set[str]]
    haplotypes: Union[List[int], List[str]]
    min_pvalue_str: str = ""
    pvalues: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    def print_members_list(self) -> str:
        """Returns a string that has all of the members ids separated by space

        Returns
        -------
        str
            returns a string where the members list attribute
            is formatted as a string for the output file. Individuals strings are joined by comma.
        """
        return ", ".join(list(map(str, self.members)))

    def __lt__(self, comp_class: T) -> bool:
        """Override the less than method so that objects can be sorted in
        ascending numeric order based on cluster id.

        Parameters
        ----------
        comp_class : Network
            Network object to compare.

        Returns
        -------
        bool
            returns true if the self cluster id is less than the
            comp_class cluster id.
        """

        return self.clst_id < comp_class.clst_id
