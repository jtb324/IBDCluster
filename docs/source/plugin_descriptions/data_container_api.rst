Structure of the Data container
===============================

There is a specific class that the clustering module of DRIVE creates. This class is called "Data" and is passed to each plugin (see code block below). This class has no methods. The following section will describe the class attributes so the user understands what information each plugin can access by default.


.. code:: python

    @dataclass
    class Data:
        """main class to hold the data from the network analysis and the different pvalues"""

        networks: List[Network_Interface]
        output_path: Path
        carriers: Dict[str, Dict[str, List[str]]]
        phenotype_descriptions: Dict[str, Dict[str, str]]

.. admonition:: Changing the Data class attributes

    The attributes described in this section are what DRIVE comes with by default. Users can add additional attributes through custom plugins. The advantages/disadvantages of customizing the class attributes will be discussed in further detail in a future section on creating custom plugins.

Data class attributes:
----------------------
- **networks**: This is a list of all the networks identified by the DRIVE clustering analysis. This list consists of Network objects that are shown by the code block below. These objects have information about the individuals in the cluster, the haplotypes in the cluster, and the size of the cluster. The Network class has attributes for pvalues, and min_pvalue_str but these values are empty until they are calculated by the "pvalues" plugin.

.. code:: python

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

.. role:: python(code)
   :language: python

- **output**: This attribute is a path object that describes where to write the drive output to. This path will have a file prefix that drive appends a suffix to when it writes to a file. If you wish to get the parent directory you can just use :python:`network_obj.output.parent`.

- **carriers**: This attribute is a dictionary of dictionaries that tells who cases, controls, and exclusions are for each phenotype. The outer key is the phenotype id. The inner dictionary has 3 keys: "cases", "controls", "excluded". The values for each of these keys are a list that contains the individual ids for each case, control, or excluded individual, respectively.

- **phenotype_descriptions**: This attribute is another dictionary of dictionaries that has a description of each phenotype. The outer key is the phenotype id. The inner key is the string phenotype and has the phenotype description as a value.