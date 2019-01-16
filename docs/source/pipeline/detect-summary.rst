Spike detection
---------------

Once you have a config file, you can detect spikes in your recording.
Any of the following commands will detect spikes:

- ``detect`` will perform spike detection/feature extraction and save results to disk.
- ``detect-sort`` will do all of the above and additionally cluster the spikes.
- ``full`` will do all of the above and additionally pull up the curation GUI.

(See the `usage section`_ for how to invoke these commands.)

To deal with large files, detection step is performed in multiple memory loading cycles.

For each recording you specify in your config file, JRCLUST will:

#. Denoise and filter samples.
#. Compute a detection threshold based on amplitude statistics (if you have not already supplied a threshold).
#. Detect peaks, i.e., points exceeding the threshold which are also genuine turning points.
#. Merge spikes detected at multiple neighboring sites, searching over a spatiotemporal window.
   Where larger spikes are detected within this limit, weaker spikes are removed.
#. Extract spatiotemporal windows around spiking events, in both raw and filtered samples.
#. Compute low-dimensional features from the resulting waveforms.

.. _`usage section`: ../usage/index
.. _`common-average referencing`: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2666412/
