.. _pipeline:

The JRCLUST spike-sorting pipeline
==================================

The JRCLUST pipeline can be broken down into the following steps:

#. config file creation (`bootstrapping`_)
#. `spike detection/feature extraction <#spike-detection>`_
#. `clustering <#spike-clustering>`_
#. `manual curation <#cluster-curation>`_ of automatic results

.. _bootstrapping:

Bootstrapping
-------------

Your config file encapsulates the choices in parameters that you make, as well as describing the relevant probe configuration.
You can create a config file by invoking ``jrc bootstrap`` or ``jrc bootstrap /path/to/metafile.meta`` from the MATLAB command window.

If you don't specify a `SpikeGLX meta file (.meta) <https://github.com/billkarsh/SpikeGLX/blob/master/Markdown/Metadata.md>`__, you will be asked to select one.
This will collect recording-specific information from your meta file and set default parameters.
(The default parameters may be inspected on the :ref:`parameters` page.)
The location of your raw recording will be inferred from your meta file, so be sure they are similarly named (this is the SpikeGLX default) and placed together in the same directory!

If you don't have a .meta file or your .meta file is missing some data, JRCLUST requests the following information:

* **Sampling rate (Hz)**: read from ``imSampRate`` (or ``niSampRate`` for NI recordings) in .meta file
* **Number of channels in file**: read from ``nSavedChans`` in .meta file
* **μV/bit**: computed from ``imAiRangeMax``, ``imAiRangeMin`` (or ``niAiRangeMax``, ``niAiRangeMin``, and ``niMNGain`` for NI recordings) in .meta file
* **Header offset (bytes)**: set to 0 for SpikeGLX recordings since no header is stored in the .bin file
* **Data type**: select from ``int16``, ``uint16``, ``single``, or ``double`` (SpikeGLX files are saved as ``int16``)

You will also be asked to confirm your config filename and path to your raw recording.
If you have multiple recordings, full paths will be separated by commas.

.. _spike-detection:

Spike detection
---------------

Once you have a config file, you can detect spikes in your recording.
Any of the following commands will **detect** spikes:

- ``detect`` will perform spike detection/feature extraction and save results to disk.
- ``detect-sort`` will do all of the above and additionally cluster the spikes.
- ``full`` is the same as ``detect-sort``, but will also pull up the curation GUI after clustering is completed.

(See the :ref:`usage section <usage>` for how to invoke these commands.)

For each recording you specify in your config file, JRCLUST will:

#. :ref:`Denoise <denoising>` and :ref:`filter <filtering>` samples.
   Also perform `common-average referencing`_ on the filtered samples.
#. Compute a :ref:`detection threshold <compute-threshold>` from the filtered samples (if you have not already supplied a threshold).
#. :ref:`Detect peaks <peak-detection>`, i.e., points exceeding the threshold which are also genuine turning points.
#. :ref:`Merge peaks <merge-peaks>` detected at multiple neighboring sites, searching over a spatiotemporal window.
   Where larger peaks are detected within this limit, weaker peaks are removed.
#. :ref:`Extract spatiotemporal windows <extract-windows>` around spiking events, in both raw and filtered samples.
#. :ref:`Compute low-dimensional features <compute-features>` from the resulting waveforms.

A detailed description of the steps in the detection process is indexed below.

.. toctree::
   :maxdepth: 3

   detect/index

Spike clustering
----------------
After the spiking events have been detected, they must be clustered by the features extracted from them.
Any of the following commands will **sort** spikes:

- ``sort`` will cluster spiking events using spikes you have detected previously.
  If JRCLUST can't find a previous detection, it will also detect them for you.
- ``detect-sort`` will cluster spiking events after detecting them.
  If you have previously detected spikes, JRCLUST will overwrite them.
- ``full`` is the same as ``detect-sort``, but will also pull up the curation GUI after clustering is completed.

(See the :ref:`usage section <usage>` for how to invoke these commands.)

JRCLUST will cluster the spike features using a variant of the clustering algorithm of `Rodriguez and Laio`_,
also known as density-peak clustering or :math:`\rho`-:math:`\delta` clustering.
The general algorithm computes pairwise distances in the feature space :math:`\mathbb{R}^n`, and, given a cutoff distance :math:`d`,
assigns to each point :math:`x_i` in the feature space a *density*

.. math::

    \rho_i := \sum_{j} I(\|x_i - x_j\|_2 < d),

where :math:`I` is the indicator function, giving 1 if the condition is true and 0 otherwise.
In plain English, :math:`\rho_i` is the number of points within a ball of radius :math:`d` centered at :math:`x_i`.

Once each point has been assigned a density, we then find the distance to the nearest neighbor of higher density,

.. math::

    \delta_i := \min_{\rho_i < \rho_j} \|x_i - x_j\|_2.

(If there is no :math:`j` such that :math:`\rho_i < \rho_j`, i.e., :math:`x_i` is a maximally dense point,
then :math:`\delta_i` is defined to be :math:`\max_j \|x_i - x_j\|_2`.)

If the point :math:`x_i` is both sufficiently dense (i.e., :math:`\rho_i` is large enough)
**and** sufficiently far from any other point with higher density, then :math:`x_i` is deemed a *cluster center*,
and all points with :math:`x_i` as nearest neighbor of higher density will be assigned to the cluster centered around :math:`x_i`.

We have left the terms *sufficiently dense* and *sufficiently far* somewhat vague, to say nothing of how the cutoff distance :math:`d` is determined.
These, along with the variations JRCLUST makes to accommodate the large volume of data, are elaborated `here <./sort/index.html>`__.

.. toctree::
   :maxdepth: 2

   sort/index

Cluster curation
----------------

Manual verification and correction of the automatic clustering can be done using the GUI curation tool.

Any of the following commands will **curate** spike clusters:

- ``manual`` will pull up the curation GUI.
  If JRCLUST can't find a previous clustering, it will ask if you want to cluster your spikes.
- ``detect-sort`` will do all of the above and additionally cluster the spikes.
- ``full`` will detect and cluster spikes, pulling up the curation GUI after clustering is completed.
  If you have previously detected spikes, JRCLUST will overwrite them.

(See the :ref:`usage section <usage>` for how to invoke these commands.)

In addition to deleting, splitting, or merging clusters, you may also annotate them as noise, MUA (multi-unit activity),
or provide arbitrary notes on clusters.
For a good primer on why you might want to do these things, take a look at
`Phy's documentation <https://phy-contrib.readthedocs.io/en/latest/template-gui/#basic-concepts>`__.
A detailed description of the different views onto the data is indexed below.

.. toctree::
   :maxdepth: 2

   curate/index

.. _`common-average referencing`: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2666412/
.. _`Rodriguez and Laio`: http://science.sciencemag.org/content/344/6191/1492
