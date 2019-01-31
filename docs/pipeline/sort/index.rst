.. _pipeline-sort:

Spike clustering
================

.. _rl-adjustments:

Adjustments to the algorithm
----------------------------

Because many millions of spikes may be detected in a single recording and the
distance matrix grows as :math:`O(n^2)`, we quickly run out of memory.
Instead, we take advantage of the fact that two spikes which are far away from each
other in space and in time are unlikely to be from the same neuron and thus should
not belong to the same cluster.
So instead of computing the distance between a spike and *every* other spike, we
restrict our comparison of spikes on site :math:`K` to spikes which either
also peak on site :math:`K`, or have a :ref:`secondary peak <compute-features>`
on site :math:`K`.
Additionally, :ref:`if you specify <nClusterIntervals>`, we also restrict our comparison
to spikes occurring nearby in time (you might want to do this if you observe
drift in your recording).
This has a tendency to oversplit clusters, but they can be
:ref:`merged later <merge-post-hoc>`.

.. _dist-cut:

Determining :math:`d`
---------------------

The best cutoff distance in feature space is generally not known *a priori*.
In order to determine a good cutoff distance, we instead estimate a
:ref:`percentile you specify <distCut>` (the 2nd percentile by default),
call it :math:`d_i`, of the pairwise distances between spike features on each
site :math:`i`.
We then take the median of the :math:`d_i` as our global cutoff distance.

(This behavior can be overridden by setting the parameter :ref:`useGlobalDistCut`
to 0 in your config file.)

.. _assign-clusters:

Assigning clusters: "sufficiently dense" and "sufficiently far"
---------------------------------------------------------------

Once every spike has its :math:`\rho` and :math:`\delta` scores, along with the
nearest neighbor of higher density corresponding to the :math:`\delta` score, we
then determine cluster centers.

You must :ref:`specify <log10RhoCut>` the :ref:`thresholds <log10DeltaCut>` to
determine a cluster center.
In particular, a spike :math:`i` must have :math:`\log_{10}(\rho_i) > \text{log10RhoCut}` ("sufficiently dense")
and :math:`\log_{10}(\delta_i) > \text{log10DeltaCut}` ("sufficiently far").
If :ref:`you specify <RDDetrendMode>`, this comparison may be made against a
detrended plot, as in the images below.

.. image:: /.static/rdPlot.png

Once cluster centers are assigned, it is a straightforward operation to assign
spikes to clusters: for each spike which has not been assigned a cluster (initially,
this is all spikes which are not cluster centers), assign it to the same cluster
as its nearest neighbor of higher density.
Most spikes will remain unassigned at first, but on each iteration of this procedure,
the number of spikes which are unassigned will shrink until all spikes have been
assigned to clusters.

.. note::
   Because we :ref:`compared primary peaks with secondary peaks <rl-adjustments>`,
   a spike may have as its nearest neighbor a spike on a nearby site.

.. _merge-post-hoc:

Merging possibly oversplit clusters
-----------------------------------

After the initial clustering using, units with similar mean raw waveforms are
merged based on their similarity, as computed using the
`Pearson correlation coefficient`_ of the waveforms.
Since probe drift changes the spike waveforms, averaged waveforms are computed
at three depth ranges per unit using their center-of-mass positions.
Additionally, a time shift is applied between each pair to find the maximum cross-correlation value.
The similarity score is taken as the maximum score from the resulting
:math:`3 \times 3` correlation matrix for each unit pair.
Any cluster pair whose maximum similarity exceeds a threshold
:ref:`you specify <maxUnitSim>` is merged.

This process is repeated :ref:`some number of times <nPassesMerge>` and the resulting clustering
is output for inspection.

.. _`Pearson correlation coefficient`: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
