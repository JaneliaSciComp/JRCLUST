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

.. Clustering step is performed using [[sort_()]] function ("sort" command) that uses a density-based clustering ([[DPCLUS]]) method (Rodriguez and Laio, Science'14).
.. [[DPCLUS]] computes two clustering parameters (density and distance) locally in time and space.
.. Distributed spatiotemporal sorting combined with GPU usage significantly reduces the clustering time.
.. [[DPCLUS]] algorithm assigns a nearest neighbor to each spike that points to a higher density gradient.
.. Cluster memberships are first assigned to the density peaks having outlier density and distance values.
.. The membership information is then recursively copied to their nearest neighbors by following the density gradient from the peak.
.. Since the nearest neighbor can exist at a different site, the membership assignment is globally propagated to all sites.
.. After the initial clustering using [[DPCLUS]], clusters with similar mean raw waveforms are merged based on the similarity score using cross-correlation.
.. To deal with probe drift, three copies of waveforms are computed per unit based on their center of mass locations.
.. The similarity score takes the highest value of the waveform pairs (3x3 matrix for each pair of units being compared).
.. The similarity score is computed by taking the highest correlation value by shifting time between two unit waveforms up to +/-0.25 ms ([[spkRefrac_ms]] parameter).

.. _`Rodriguez and Laio`: http://science.sciencemag.org/content/344/6191/1492
