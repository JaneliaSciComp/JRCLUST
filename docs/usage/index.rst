.. _usage:

Using JRCLUST
=============

.. toctree::
   :glob:
   :maxdepth: 2
   :caption: Table of contents
   :name: usagetoc

   tutorial
   migrating
   io

Getting started
---------------

You've got a recording file and you've got MATLAB.
Now what?
If you've used an older version of JRCLUST, you might want to check out the :ref:`migrating <migrating>` page to see what's changed from the old JRCLUST.
If you're completely new to JRCLUST, start with the :ref:`tutorial <tutorial>`.

System requirements
-------------------

Software Requirements
~~~~~~~~~~~~~~~~~~~~~

JRCLUST 4 was built and tested with **MATLAB R2017b** on Microsoft **Windows 10**.
If you are running another operating system or an earlier version of MATLAB, JRCLUST may not work as expected.
Post an issue on our `issues page`_ with the error message and we'll get it fixed.

You will also need the `CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit-archive>`__.
See `this link <https://www.mathworks.com/help/distcomp/gpu-support-by-release.html>`__ to determine which CUDA Toolkit version you will need for your version of MATLAB.
It's best to keep your GPU driver updated as well.

In addition, the following toolboxes are **required**:

- `Parallel Computing <https://www.mathworks.com/products/parallel-computing.html>`__
- `Image Processing <https://www.mathworks.com/products/image.html>`__
- `Signal Processing <https://www.mathworks.com/products/signal.html>`__
- `Statistics and Machine Learning <https://www.mathworks.com/products/statistics.html>`__

Hardware Requirements
~~~~~~~~~~~~~~~~~~~~~

An NVIDA GPU with `compute capability`_ 3.5 or later is strongly recommended.
In addition, we also recommend you have at least 32 GiB of memory.

.. _usage-cli:

Usage
-----

In general, usage goes like

.. code-block:: matlab

   jrc COMMAND ARG1 ARG2 ...

where ``COMMAND`` is described below.

Documentation and help
~~~~~~~~~~~~~~~~~~~~~~

-  ``jrc help``: Display a help menu.
-  ``jrc version``: Display the version number.
-  ``jrc about``: Display program information.

Pipeline commands
~~~~~~~~~~~~~~~~~

-  ``jrc importv3 mysession_jrc.mat``: Import a JRCv3 ``_jrc.mat`` file to the new format.
   **Does not overwrite your old results.**
-  ``jrc bootstrap``: Create a :ref:`config file <bootstrapping>`.
-  ``jrc detect myparams.prm``: Perform spike detection and feature extraction.
   (Replace ``myparams`` with a unique name of your choosing.)
   Output files:

   - myparams\_res.mat: A MAT-file containing detection results
   - myparams\_filt.jrc (filtered spike traces)
   - myparams\_raw.jrc (raw spike traces)
   - myparams\_features.jrc (computed features)

   See :ref:`io-files` for details.

-  ``jrc sort myparams.prm``: Cluster detected spike features.
   If JRCLUST can't find any detected spikes, it will also perform the detect step.
   Your myparams\_res.mat file will be updated with a :ref:`clustering handle <DensityPeakClustering>` object.
-  ``jrc detect-sort myparams.prm``: Detect and cluster spikes.
   The difference between ``sort`` and ``detect-sort`` is that ``detect-sort`` will overwrite any previously-detected spikes,
   whereas ``sort`` will only perform spike detection if JRCLUST can't find any previously-detected spikes.
-  ``jrc recluster myparams.prm``: Recluster spikes.
   You might want to do this if, e.g., you've updated the clustering threshold parameters but don't want to recompute rho and delta.
-  ``jrc manual myparams.prm``: Open the manual curation GUI.

Preview commands
~~~~~~~~~~~~~~~~

-  ``jrc probe {myprobe.prb, myparams.prm}``: Plot probe layout.
   You may specify either a config file or a probe file directly.
-  ``jrc traces myparams.prm [File #]``: Display traces from your raw recordings.
   If you have multiple recordings, you may specify the index of the file to display from.
-  ``jrc activity myparams.prm``: Plot firing rate as a function of time and depth.

.. _missing-commands:

Missing commands
~~~~~~~~~~~~~~~~

Many other commands were available with previous versions of JRCLUST.
These have been disabled because, as far as we know, they were not widely used.
If you're missing a command, post an issue on our `issues page`_ request and we
will correct our error.

Getting help
------------

Have questions?
Check out our support channel on `Gitter <https://gitter.im/JRCLUST/community>`__.
If you run into any bugs or have any feature requests, `create an issue <https://github.com/JaneliaSciComp/JRCLUST/issues/new>`__ on our GitHub `issues page`_.

.. _`issues page`: https://github.com/JaneliaSciComp/JRCLUST/issues
.. _`compute capability`: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compute-capability
