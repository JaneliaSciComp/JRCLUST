.. _parameters:

JRCLUST parameters
==================

Common parameters
-----------------

.. _CARMode:

``CARMode``
^^^^^^^^^^^

The meaning of "average" in "common average reference".

One of the following:

- ``mean``: (**default**) Computes the mean across sites, excluding sites in :ref:`ignoreSites`.
- ``median``: Computes the median across sites, excluding sites in :ref:`ignoreSites`.
- ``none``: Skips CAR.

.. _clusterFeature:

``clusterFeature``
^^^^^^^^^^^^^^^^^^

The feature to extract from your spike waveforms in order to cluster them.

One of the following:

- ``gpca``: :ref:`"Global" PCA <gpca>`.
   :ref:`nPcPerChan` sets the number of components per channel (default: 1).
   .. warning:: This feature is currently disabled due to an issue in formulation.
- ``pca``: Projection of spike waveforms onto :ref:`principal components <pca>`.
- ``cov``: Covariance feature.
- ``vpp``: Peak-to-peak spike amplitudes.
- ``vmin``: Minimum spike amplitude
- ``vminmax``: Minimum and maximum spike amplitudes (2 features/site)
- ``energy``: Standard deviation of spike waveform

.. _evtDetectRad:

``evtDetectRad``
^^^^^^^^^^^^^^^^

Maximum distance (in μm) for :ref:`extracting spike waveforms <extract-windows>` (``r2`` in the figure below).

**Default**: 50

.. image:: /.static/evtDetectRad.png
   :scale: 25%

.. _evtMergeRad:

``evtMergeRad``
^^^^^^^^^^^^^^^

Maximum distance (in μm) to search over for :ref:`potential duplicates <merge-peaks>` (``r1`` in the figure below).
This distance is used to determine the number of sites to extract features if :ref:`nSiteDir` is empty.

**Default**: 75

.. image:: /.static/evtDetectRad.png
   :scale: 25%

.. _evtWindow:

``evtWindow``
^^^^^^^^^^^^^

Time range (in ms) of filtered spike waveforms, centered at the peak.

Must be an array with 2 elements, the first negative and the second positive.
For example, if ``evtWindow`` is set to [-0.5, 0.5], then 1/2 ms worth of samples
are extracted before and after the spiking event.

**Default**: [-0.25, 0.75]

.. _evtWindowRaw:

``evtWindowRaw``
^^^^^^^^^^^^^^^^

Time range (in ms) of raw spike waveforms, centered at the peak.

Must be an array with 2 elements, the first negative and the second positive.
For example, if ``evtWindowRaw`` is set to [-1, 1], then 1 ms worth of samples
are extracted before and after the spiking event.

**Default**: [-0.5, 1.5]

.. _fftThresh:

``fftThresh``
^^^^^^^^^^^^^^^^

The threshold, in units of `MAD`_, to apply to the power-frequency product of your samples when :ref:`denoising <denoising>`.
Frequencies with power-frequency product above this threshold will be zeroed out as noise.

Must be a nonnegative number.
The recommended value is 10.
Setting to 0 disables the notch filter (**default** behavior).

.. _filterType:

``filterType``
^^^^^^^^^^^^^^

The type of filter to apply to your raw samples.

One of the following:

- ``ndiff``: (**default**) Applies a differentiation filter a kernel depending on the order given in :ref:`nDiff_filt`.
- ``sgdiff``: Applies a `Savitzky-Golay <https://en.wikipedia.org/wiki/Savitzky–Golay_filter>`_ filter depending on the order given in :ref:`nDiff_filt`.
- ``bandpass``:
- ``fir1``:
- ``user``: Convolves your raw samples with a kernel of your choosing.
- ``none``: Skips filtering.

.. _ignoreSites:

``ignoreSites``
^^^^^^^^^^^^^^^

.. _log10DeltaCut:

``log10DeltaCut``
^^^^^^^^^^^^^^^^^

.. _log10RhoCut:

``log10RhoCut``
^^^^^^^^^^^^^^^

.. _maxWavCor:

``maxWavCor``
^^^^^^^^^^^^^

.. _nDiff_filt:

``nDiff_filt``
^^^^^^^^^^^^^^

.. _nSiteDir:

``nSiteDir``
^^^^^^^^^^^^

The number of neighboring sites to group in either direction.
The total number of sites per spike group (``nSitesEvt``) is 1 + 2\*``nSiteDir``.
In other words, a spike group includes the site on which the spike occurs, along with ``nSiteDir``
sites in the horizontal direction and ``nSiteDir`` in the vertical direction.

Must be empty or a positive integer.
If empty, the number of sites per spike group is determined from :ref:`evtDetectRad`.

.. warning::
   This parameter will be deprecated in an upcoming release in favor of ``evtDetectRad``.

.. _qqFactor:

``qqFactor``
^^^^^^^^^^^^

Multiplier of the :ref:`estimate <compute-threshold>` :math:`\sigma_{\text{noise}}^{(i)}`
of standard deviation of noise distribution on each site to compute the threshold for that site.
In other words,

.. math::

    \text{Thr}_i := \text{qqFactor} \cdot \sigma_{\text{noise}}^{(i)}

is the spike detection threshold for site :math:`i`.

Must be a positive number.
**Default** is 5.

Advanced parameters
-------------------

The following parameters can be safely ignored.

.. _distCut:

``distCut``
^^^^^^^^^^^^^^^

The percentile to use when

.. _fDetectBipolar:

``fDetectBipolar``
^^^^^^^^^^^^^^^^^^

Flag to detect peaks with positive amplitudes as well as negative if true.

Must be 0 or 1.
**Default** is 0.

.. _fRealign_spk:

``fRealign_spk``
^^^^^^^^^^^^^^^^

.. _nFet_use:

``nFet_use``
^^^^^^^^^^^^

.. _nPassesMerge:

``nPassesMerge``
^^^^^^^^^^^^^^^^

.. _nneigh_min_detect:

``nneigh_min_detect``
^^^^^^^^^^^^^^^^^^^^^

Number of sample neighbors exceeding the detection threshold required for a sample
to be considered a peak.
For example, consider a putative peak occurring at sample :math:`t_i`.
If ``nneigh_min_detect`` is set to 1, then **either** sample :math:`t_{i-1}` or :math:`t_{i+1}`
must also exceed the detection threshold.
If ``nneigh_min_detect`` is set to 2, then **both** sample :math:`t_{i-1}` and :math:`t_{i+1}`
must also exceed the detection threshold.
If ``nneigh_min_detect`` is set to 0, then samples :math:`t_{i-1}` and :math:`t_{i+1}`
are not considered.

Must be one of 0, 1, or 2.
**Default** is 0.

.. _nPcPerChan:

``nPcPerChan``
^^^^^^^^^^^^^^

Number of principal components to take per site.
**Default** is 1.

.. _nTime_clu:

``nTime_clu``
^^^^^^^^^^^^^

When clustering, take the :math:`\frac{1}{\text{nTime_clu}}` fraction of all
spikes around a spiking event to compute distance.

For example, if ``nTime_clu`` = 1, all spikes will be used;
if ``nTime_clu`` = 2, JRCLUST will take the half of all spikes which are closest
in time to compute distances.
Increasing this value will take fewer and fewer spikes to compare at the risk of
oversplitting clusters (you might want to do this if you observe fast drift in your
recording).
However, automated merging based on the :ref:`waveform correlation <maxWavCor>`
can merge most of the units initially split by drift.

**Default** is 4.

.. _refracIntms:

``refracIntms``
^^^^^^^^^^^^^^^

Spike refractory period, in milliseconds.
**Default** is 25.

.. _rlDetrendMode:

``rlDetrendMode``
^^^^^^^^^^^^^^^^^

.. _useGlobalDistCut:

``useGlobalDistCut``
^^^^^^^^^^^^^^^^^^^^

Flag to estimate a global :ref:`distance cutoff <dist-cut>` in the clustering algorithm.
Set to 0 to use a separate distance cutoff for each site.

**Default** is 1.

.. _userFiltKernel:

``userFiltKernel``
^^^^^^^^^^^^^^^^^^

Custom filter kernel.
Your filtered samples will be the output of a convolution of your raw samples
with this kernel.
You must specify this if and only if your :ref:`filterType` is ``'user'``.

**Default** is empty.

.. _MAD: https://en.wikipedia.org/wiki/Median_absolute_deviation
