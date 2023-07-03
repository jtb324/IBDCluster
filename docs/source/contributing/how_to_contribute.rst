How to contribute
=================

Aliqua elit cillum irure nulla commodo occaecat nisi adipisicing non magna tempor. Adipisicing ex minim dolor sint reprehenderit ex sit veniam id excepteur do voluptate. Excepteur incididunt nisi consequat laboris. Cillum irure sunt sit culpa anim occaecat do cupidatat culpa deserunt. Esse adipisicing non esse amet eiusmod. Laborum id incididunt ut Lorem ea ea consectetur ut exercitation ullamco nulla laboris. Et laborum reprehenderit in ullamco ea ullamco duis et.

Improving Documentation:
------------------------
One of the key challenges in the scientific community is documentation. Oftentimes there is a lack of updated documentation, there is a lack of any documentation, or the documentation is not suitable for individuals that have not used the software before. Scientific documentation often lacks the technical detail required to run the program. One of the most helpful contributions can be to either propose a change to the documentation that broadly falls into one of the following four categories:

    1. Add missing information that is needed to successfully run the program and get accurate results.
    2. Expand upon information in the documentation. (Information in the documentation may be incomplete because of scientific advances) 
    3. Update information that may no longer be correct due to scientific advances. 
    4. Improve clarity of the documentation

Improving Code:
---------------
Due to changes in the python ecosystem the code in DRIVE will need to be updated to prevent technical debt from piling up. One of DRIVE's goal is to maintain backwards compatibility but it also has to remain usable meaning development of the tool has to walk a tight line between implementing new features and also not making the tool harder for the general scientific community to use. Contributors are expected to run unit test to ensure that the proposed changes produce the same results. In cases where it is deemed appropriate to break the current DRIVE api, contributors are expected to use `semantic versioning <https://semver.org/>`_ to indicate that the change is API breaking. Contributors are also expected to format code using pylint, black, and isort to maintain code style consistency. 