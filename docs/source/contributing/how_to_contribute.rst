How to contribute
=================

Fixing Bugs:
------------
If you discover any bugs, you could either submit and issue through Github using the recommended issue template. You are also welcome to try to fix the issue on your own and then submit a pull request. We ask that you read the code style section (:ref:`code_style`). Your pull request will only be accepted if it conforms to the development standards that DRIVE uses.

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