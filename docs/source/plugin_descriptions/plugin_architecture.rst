How is DRIVE extensible
=======================

DRIVE is designed with user extensibility in mind and it accomplishes this through a plugin architecture. This design allows user to dynamic load their own code at runtime (as long as the code conforms to the specified Plugin Interface). You can read more about this design strategy at this site `Plugin Architecture <dotcms.com/blog/post/plugin-achitecture>`_. The only component of DRIVE that can be changed is the clustering algorithm. Everything else, such as the statistics and how the program writes to an output file are completely customizable to the user.

At its core, DRIVE, is just a clustering software and then users can bring their own extensions to customize it to fit their needs. That being said DRIVE, offers two plugins out of the box. One plugin uses a binomial test to calculate phenotypic enrichment within each network and the other plugin writes the network information to a file. This plugins are described in more detail in this section:

- :doc:`factory_plugins`
