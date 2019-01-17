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
#. Compute low-dimensional features from the resulting waveforms.

JRCLUST will chunk up large files and perform all these steps in sequence over each chunk.

.. _`common-average referencing`: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2666412/
