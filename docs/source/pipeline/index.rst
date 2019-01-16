The JRCLUST spike-sorting pipeline
==================================

The JRCLUST pipeline can be broken down into the following steps:

#. config file creation (bootstrapping)
#. spike detection/feature extraction
#. clustering
#. manual curation of automatic results

.. toctree::
  :maxdepth: 1
  :hidden:

  bootstrap-summary
  detect-summary
  detect/index

.. _bootstrap:
.. include:: ./bootstrap-summary.rst

.. _detect:
.. include:: ./detect-summary.rst
..
.. Clustering step is performed using [[sort_()]] function ("sort" command) that uses a density-based clustering ([[DPCLUS]]) method (Rodriguez and Laio, Science'14). [[DPCLUS]] computes two clustering parameters (density and distance) locally in time and space. Distributed spatiotemporal sorting combined with GPU usage significantly reduces the clustering time. [[DPCLUS]] algorithm assigns a nearest neighbor to each spike that points to a higher density gradient. Cluster memberships are first assigned to the density peaks having outlier density and distance values. The membership information is then recursively copied to their nearest neighbors by following the density gradient from the peak. Since the nearest neighbor can exist at a different site, the membership assignment is globally propagated to all sites. After the initial clustering using [[DPCLUS]], clusters with similar mean raw waveforms are merged based on the similarity score using cross-correlation. To deal with probe drift, three copies of waveforms are computed per unit based on their center of mass locations. The similarity score takes the highest value of the waveform pairs (3x3 matrix for each pair of units being compared). The similarity score is computed by taking the highest correlation value by shifting time between two unit waveforms up to +/-0.25 ms ([[spkRefrac_ms]] parameter).
..
.. Manual verification and correction can be carried out using GUI tool ("manual" command, see [[manual_()]]). The GUI tool is comprised of nine linked views that are interactively updated by user actions.
..
.. # [[detect_()]]: filter, detect spikes, extract features
.. * The process is performed in multiple memory loading cycles to deal with large files. Each memory loading step calls [[file2spk_()]] by passing the raw traces loaded from the recording file.
.. * Output files: filtered spike waveforms are saved to [[_spkwav.jrc]], raw spike waveforms are saved to [[_spkraw.jrc]], and clustering features are saved to [[_spkfet.jrc]]. The rest of the files are saved to [[_jrc.mat]] in matlab file format. Each output file uses the parameter file names as a prefix (excluding ".prm" extension).
.. * [[file2spk_()]]: load recordings from file(s) and divide the problem to the available memory size and merge the solutions
..     * [[wav2spk_()]]: Filter, detect and group spikes, and extract features for each memory loading block
..         * [[fft_clean()]]: Adaptive notch filter (skip if P.fft_thresh==0)
..         * [[filt_car_()]]: Filter (P.[[vcFilter]]) and apply common average referencing (P.vcCommonRef)
..         * [[detect_spikes_()]]: Detect spikes using a threshold setting (P.[[qqFactor]]).
..             * [[spikeDetectSingle_fast_()]]: Perform spike detection for each channel using GPU
..             * [[spikeMerge_()]]: Merge spikes detected at multiple neighboring sites using multiple CPU cores
..         * [[mn2tn_wav_()]]: Extract spike waveforms from a list of specified site (viSite_spk) and time (viTime_spk)
..         * [[trWav2fet_()]]: Extract spike features from the supplied spike waveforms (see [[P]].[[vcFet]])
..
.. ![""](https://github.com/JaneliaSciComp/JRCLUST/blob/master/img/event_group.png)
..
.. # [[sort_()]]: cluster features and merge based on waveform similarities
.. * sort spikes using the features ([[fet2clu_()]]), and automatically merge similar clusters using waveform similarities ([[post_merge_()]]).
.. * Output of the clustering is stored as "S_clu" struct in the master struct ([[S0]]), which is saved to ([[_jrc.mat]]).
.. * [[fet2clu_()]]: Cluster the spike features ([[trFet_spk]]) using [[DPCLUS]] clustering algorithm (Rodriguez-Laio, Science'14).
..     * [[cluster_spacetime_()]]: Apply spatiotemporal [[DPCLUS]] (divide and conquer density-based clustering)
..         * For each site, compute Rho (density) and Delta (distance) values for each spike using [[DPCLUS]] clustering
..         * [[calc_dc2_()]]: Compute distance cut-off value (dc) based on the P.[[dc_percent]] parameter.
..         * [[cuda_rho_()]]: Compute local density (rho) per each spike using GPU.
..         * [[cuda_delta_()]]: Compute separation distance (delta) per each spike using GPU.
..     * [[postCluster_()]]: Group clusters using rho and delta values (see P.[[rho_cut]], P.[[delta1_cut]]), and apply automated merging based on the mean unit waveforms
..     * [[post_merge_()]]: Merge units based on the mean spike waveforms. Unit pairs with correlation above the threshold (P.maxWavCor) is automatically merged.
..         * [[post_merge_wav_()]]: Apply automated merging based on the waveform similarity
..         * [[S_clu_position_()]]: Compute cluster positions by computing the center of mass
..         * [[S_clu_quality_()]]: Compute the cluster quality scores
..
.. ![""](https://github.com/JaneliaSciComp/JRCLUST/blob/master/img/rho_delta_plot.png)
..
.. # [[manual_()]]: manual spike sorting GUI
.. * Manual verification and correction can be carried out using GUI tool ("manual" command, see [[manual_()]]). The GUI tool is comprised of nine linked views that are interactively updated by user actions.
.. * [[plot_FigPos_()]]: Plot or update the position view (spike waveforms arranged by the site positions)
.. * [[plot_FigMap_()]]: Plot or update the color-coded activity map on the probe site layout
.. * [[plot_FigWav_()]]: Plot or update the average unit waveform view
.. * [[plot_FigTime_()]]: Plot or update the amplitudes vs. time view
.. * [[plot_FigProj_()]]: Plot or update the 2D projection view
.. * [[plot_FigWavCor_()]]: Plot or update the waveform correlation matrix
.. * [[plot_FigHist_()]]: Plot or update the ISI histogram
.. * [[plot_FigIsi_()]]: Plot or update the [ISI return map](https://en.wikipedia.org/wiki/Poincar%C3%A9_map) (previous ISI vs. next ISI).
.. * [[plot_FigCorr_()]]: Plot or update the spike timing cross-correlation view
.. * [[plot_FigRD_()]]: Plot the Rho-Delta plot (Rodriguez-Laio) and identify the cluster centers that exceeds the cluster detection threshold
