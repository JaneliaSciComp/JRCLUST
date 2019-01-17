JRCLUST parameters
==================

Common parameters
-----------------

.. _carMode:

``carMode``
^^^^^^^^^^^

The meaning of "average" in "common average reference".

One of the following:

- ``mean``: (**default**) Computes the mean across sites, excluding sites in :ref:`ignoreSites`.
- ``median``: Computes the median across sites, excluding sites in :ref:`ignoreSites`.
- ``none``: Skips CAR.

.. _evtDetectRad:

``evtDetectRad``
^^^^^^^^^^^^^^^^



.. _fftThreshMAD:

``fftThreshMAD``
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
- ``fir1``
- ``user``
- ``none``: Skips filtering.

.. _ignoreSites:

``ignoreSites``
^^^^^^^^^^^^^^^

.. _nDiff_filt:

``nDiff_filt``
^^^^^^^^^^^^^^

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

.. _fDetectBipolar:

``fDetectBipolar``
^^^^^^^^^^^^^^^^^^

Flag to detect peaks with positive amplitudes as well as negative if true.

Must be 0 or 1.
**Default** is 0.

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

.. _refracIntms:

``refracIntms``
^^^^^^^^^^^^^^^

Spike refractory period, in milliseconds.
**Default** is 25.

.. _MAD: https://en.wikipedia.org/wiki/Median_absolute_deviation
