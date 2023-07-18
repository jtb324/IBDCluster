Testing
=======

DRIVE uses the `Pytest framework <https://docs.pytest.org/en/7.3.x/>`_ to ensure that code is maintaining backward's compatibility and performing as expected. This process is done throught unit test or integration tests. These test are automatically run using Github actions when a proposed change is pushed to Github so there is no reason for a contributor to run these test. That being said, if you wish to run the tests you can use either of the following commands:

.. code::

    pytest -v -m unit

    or 

    pytest -v -m integtest

Adding Tests:
------------

DRIVE is not fully covered by unit test so a very helpful way to contribute to drive would be to add unit test. All we ask is that you use the pytest mark "unit" if you are adding a unit test or "integtest" if you are adding an integration test. The following section of the Pytest documentation can explain this process in more detail: `How to create tests <https://docs.pytest.org/en/7.3.x/getting-started.html>`_ .


