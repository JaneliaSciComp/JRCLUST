.. _ClusteringInterface:

The Clustering interface
------------------------

This section is mainly about extending the JRCLUST framework to
handle different clustering algorithms or the results of other
spike sorting packages.
The first thing to do is to design a class that implements the ``Clustering``
interface.

The ``Clustering`` interface describes a clustering of spikes.
It can be found in the ``jrclust.interfaces.Clustering`` class
(i.e., in +jrclust/@Clustering/Clustering.m).

.. _Clustering-props:

Clustering properties
~~~~~~~~~~~~~~~~~~~~~

hCfg
++++

A :ref:`Config` object.

.. _DensityPeakClustering:

DensityPeakClustering
---------------------

Coming soon.
