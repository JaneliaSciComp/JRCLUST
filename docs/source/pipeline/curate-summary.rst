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
