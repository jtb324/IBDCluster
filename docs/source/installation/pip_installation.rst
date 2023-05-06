Installing DRIVE from Pip
=========================
DRIVE v1.0.0 has been officially released on `PYPI <https://pypi.org/project/drive-ibd/>`_! This installation method is recommended for those who wish to use the software without any modification. PIP will install all of the necessary dependencies of DRIVE so that the user doesn't have to worry about dependency management. 

.. note::

    The recommend way to install DRIVE using pip would be to either create a virtual environment using `Anaconda <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ or `venv <https://docs.python.org/3/library/venv.html>`_. Once you activate the virtual environment then you can use the above pip command to install DRIVE into an isolated environment.


DRIVE can be installed using the following command:

.. code::

    pip install drive-ibd

.. warning::
    
    DRIVE has been tested with python ranging from v3.6.13 to v3.9. Officially, DRIVE v1.0.0 only supports python 3.8 and 3.9 because the documentation tool sphinx only supports python >= 3.8. If the user only installs the basic dependencies (meaning not the docs dependency group) then they can use python 3.6 or 3.7 if so desired.
    