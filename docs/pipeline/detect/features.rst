.. _feature-types:

Features for clustering
~~~~~~~~~~~~~~~~~~~~~~~

The following features are available for clustering your spikes.

.. _pca:

Principal components analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Principal components analysis <https://en.wikipedia.org/wiki/Principal_component_analysis>`__,
or PCA, is a method for reducing the individual observations (in our case,
a spatiotemporal window around a spiking event) to just a few components.
PCA constructs an orthogonal basis in such a way that the first basis vector explains as
much of the variance in the data as possible, and each subsequent basis vector explains less
variance, but more than any of the following basis vectors.

The diagram below illustrates the general idea:

.. image:: https://upload.wikimedia.org/wikipedia/commons/thumb/f/f5/GaussianScatterPCA.svg/1280px-GaussianScatterPCA.svg.png

Shown is a sample from a bivariate Gaussian distribution with the first two principal
components superimposed on the data points (courtesy of Wikipedia).
The larger vector describes the direction of most variance, while the smaller (orthogonal)
vector describes less than the first.
If we were to project each point onto the first principal component, we would have reduced the
dimensionality of the points in this distribution from 2 to 1.
This is a trivial example, but if we extend the 2-dimensional Gaussian to :math:`n` dimensions
we see how useful PCA can be in dimensionality reduction.

If you select 'pca' as your :ref:`clusterFeature`, then
JRCLUST computes principal components from a random subsample of at most 10,000 of the :ref:`traces <extract-windows>`
of detected spikes in a given :ref:`chunk of data <chunking>`.
All spike traces on each site in the event radius are then projected onto the first
(:ref:`or first 2 or 3 <nPCsPerSite>`) principal vectors, giving us 1 (or 2, or 3) points per site per spike,
instead of the entire waveform.

Additionally, :ref:`if you specify <interpPC>`, JRCLUST will interpolate the first principal component in order
to maximize the projection of the spike traces onto that component.

.. _gpca:

Global PCA
^^^^^^^^^^

Global PCA describes the same procedure as PCA, :ref:`above <pca>`, but instead using a single common set of
principal vectors computed from the first chunk, rather than computing a set of principal vectors
for each chunk.
You can use global PCA as your feature by setting your :ref:`clusterFeature` to 'gpca'.

.. _vmin:

Vmin
^^^^

This feature projects the traces in the event window of each spike onto their minimum values.

.. _vminmax:

Vminmax
^^^^^^^

This feature projects the traces in the event window of each spike onto their minimum **and** maximum values.
This gives 2 features per site per spike.

.. _feature-vpp:

Vpp
^^^

This feature projects the traces in the event window of each spike onto their minimum and maximum values,
taking the distance between them as the feature.
