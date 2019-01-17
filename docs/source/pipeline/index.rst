The JRCLUST spike-sorting pipeline
==================================

The JRCLUST pipeline can be broken down into the following steps:

#. config file creation (`bootstrapping`_)
#. `spike detection/feature extraction <#spike-detection>`_
#. `clustering <#spike-clustering>`_
#. `manual curation <#cluster-curation>`_ of automatic results

Bootstrapping
-------------
.. include:: ./bootstrap-summary.rst

Spike detection
---------------
.. include:: ./detect-summary.rst
A detailed description of the steps in the detection process is `here <./detect/index.html>`__.

Spike clustering
----------------
.. include:: ./sort-summary.rst
These, along with the variations JRCLUST makes to accommodate the large volume of data, are elaborated `here <./sort/index.html>`__.

Cluster curation
----------------
.. include:: ./curate-summary.rst
A detailed description of the different views onto the data is `here <./curate/index.html>`_.
