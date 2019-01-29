.. _parameters:

JRCLUST parameters
==================

Common parameters
------------------

.. _CARMode:

``CARMode``
^^^^^^^^^^^

(Formerly ``vcCommonRef``)

The meaning of 'average' in 'common average reference'.

One of the following:

- 'mean'
- 'median'
- 'none'

**Default** is 'mean'.

.. _RDDetrendMode:

``RDDetrendMode``
^^^^^^^^^^^^^^^^^

(Formerly ``vcDetrend_postclu``)

Detrending mode to apply to rho-delta values in order to determine cluster centers.

One of the following:

- 'global'
- 'local'
- 'logz'
- 'hidehiko'
- 'none'

**Default** is 'global'.

.. _autoMergeBy:

``autoMergeBy``
^^^^^^^^^^^^^^^

(Formerly ``autoMergeCriterion``)

Metric to use for automerging clusters.

One of the following:

- 'pearson'
- 'dist'

**Default** is 'pearson'.

.. _bitScaling:

``bitScaling``
^^^^^^^^^^^^^^

(Formerly ``uV_per_bit``)

ADC bit scaling factor (Conversion factor for ADC bit values to μV).

**Default** is 0.30518.

.. _blankPeriod:

``blankPeriod``
^^^^^^^^^^^^^^^

(Formerly ``blank_period_ms``)

Duration of blanking period (in ms) when the common mean exceeds blankThresh.

**Default** is 5.

.. _blankThresh:

``blankThresh``
^^^^^^^^^^^^^^^

(Formerly ``blank_thresh``)

Threshold (in MADs) above which to reject samples exceeding channel median after filtering.

**Default** is empty.

.. _clusterFeature:

``clusterFeature``
^^^^^^^^^^^^^^^^^^

(Formerly ``vcFet``)

The feature to extract from your spike waveforms in order to cluster them.

One of the following:

- 'cov'
- 'energy'
- 'pca'
- 'vmin'
- 'vminmax'
- 'vpp'

**Default** is 'pca'.

.. _dataType:

``dataType``
^^^^^^^^^^^^

(Formerly ``vcDataType``)

Format of raw recordings.

One of the following:

- 'int16'
- 'uint16'
- 'int32'
- 'uint32'
- 'single'
- 'double'

**Default** is 'int16'.

.. _dispTimeLimits:

``dispTimeLimits``
^^^^^^^^^^^^^^^^^^

(Formerly ``tlim``)

Time range (in ms) to display.

**Default** is [0, 0.2].

.. _distCut:

``distCut``
^^^^^^^^^^^

(Formerly ``dc_percent``)

Percentile of pairwise distances between spikes on a site to use as a cutoff distance.

**Default** is 2.

.. _evtDetectRad:

``evtDetectRad``
^^^^^^^^^^^^^^^^

(Formerly ``maxDist_site_spk_um``)

Maximum distance (in μm) for extracting spike waveforms.

**Default** is 75.

.. _evtWindow:

``evtWindow``
^^^^^^^^^^^^^

(Formerly ``spkLim_ms``)

Time range (in ms) of filtered spike waveforms, centered at the peak.

**Default** is [-0.25, 0.75].

.. _filtOrder:

``filtOrder``
^^^^^^^^^^^^^

Bandpass filter order.

**Default** is 3.

.. _filterType:

``filterType``
^^^^^^^^^^^^^^

(Formerly ``vcFilter``)

Type of filter to use on raw data.

One of the following:

- 'ndiff'
- 'sgdiff'
- 'bandpass'
- 'fir1'
- 'user'
- 'none'

**Default** is 'ndiff'.

.. _freqLimBP:

``freqLimBP``
^^^^^^^^^^^^^

(Formerly ``freqLim``)

Frequency cutoffs for bandpass filter.

**Default** is [300, 3000].

.. _headerOffset:

``headerOffset``
^^^^^^^^^^^^^^^^

(Formerly ``header_offset``)

Recording file header offset (in bytes).

**Default** is 0.

.. _ignoreSites:

``ignoreSites``
^^^^^^^^^^^^^^^

(Formerly ``viSiteZero``)

Sites to ignore manually.

**Default** is empty.

.. _log10DeltaCut:

``log10DeltaCut``
^^^^^^^^^^^^^^^^^

(Formerly ``delta1_cut``)

Log10 of delta cutoff (Spikes with delta values below this cutoff will not be considered as cluster centers).

**Default** is 0.6.

.. _log10RhoCut:

``log10RhoCut``
^^^^^^^^^^^^^^^

(Formerly ``rho_cut``)

Log10 of rho cutoff (Spikes with rho values below this cutoff will not be considered as cluster centers).

**Default** is -2.5.

.. _maxUnitSim:

``maxUnitSim``
^^^^^^^^^^^^^^

(Formerly ``maxWavCor``)

Threshold for merging two units having similar spike waveforms (Units with a similiarity score above this value will be merged).

**Default** is 0.98.

.. _minClusterSize:

``minClusterSize``
^^^^^^^^^^^^^^^^^^

(Formerly ``min_count``)

Minimum number of spikes per cluster (Automatically set to the maximum of this value and twice the number of features).

**Default** is 30.

.. _nChans:

``nChans``
^^^^^^^^^^

Number of channels stored in recording file (Distinct from the number of AP sites).

**Default** is 384.

.. _nClusterIntervals:

``nClusterIntervals``
^^^^^^^^^^^^^^^^^^^^^

(Formerly ``nTime_clu``)

Number of intervals to divide the recording into around a spike (When clustering, take the 1/nClusterIntervals fraction of all spikes around a spiking event to compute distance).

**Default** is 4.

.. _nPCsPerSite:

``nPCsPerSite``
^^^^^^^^^^^^^^^

(Formerly ``nPcPerChan``)

Number of principal components to compute per site.

**Default** is 1.

.. _nSiteDir:

``nSiteDir``
^^^^^^^^^^^^

(Formerly ``maxSite``)

Number of neighboring sites to group in either direction (nSitesEvt is set to 1 + 2*nSiteDir - nSitesExcl).

**Default** is empty.

.. _nSitesExcl:

``nSitesExcl``
^^^^^^^^^^^^^^

(Formerly ``nSites_ref``)

Number of sites to exclude from the spike waveform group.

**Default** is empty.

.. _nSpikesFigProj:

``nSpikesFigProj``
^^^^^^^^^^^^^^^^^^

(Formerly ``nShow_proj``)

Maximum number of spikes per cluster to display in the feature projection view.

**Default** is 500.

.. _nSpikesFigWav:

``nSpikesFigWav``
^^^^^^^^^^^^^^^^^

(Formerly ``nSpk_show``)

Maximum number of spikes per cluster to display generally.

**Default** is 30.

.. _outputDir:

``outputDir``
^^^^^^^^^^^^^

Directory in which to place output files (Will output to the same directory as this file if empty).

**Default** is an empty string.

.. _probePad:

``probePad``
^^^^^^^^^^^^

(Formerly ``vrSiteHW``)

Recording contact pad size (in μm) (Height x width).

**Default** is empty.

.. _psthTimeLimits:

``psthTimeLimits``
^^^^^^^^^^^^^^^^^^

(Formerly ``tlim_psth``)

Time range (in s) over which to display PSTH.

**Default** is empty.

.. _qqFactor:

``qqFactor``
^^^^^^^^^^^^

Spike detection threshold (Thr = qqFactor*med(abs(x-med(x)))/0.6745).

**Default** is 5.

.. _rawRecordings:

``rawRecordings``
^^^^^^^^^^^^^^^^^

Path or paths to raw recordings to sort.

**Default** is [""].

.. _refracInt:

``refracInt``
^^^^^^^^^^^^^

(Formerly ``spkRefrac_ms``)

Spike refractory period (in ms).

**Default** is 0.25.

.. _sampleRate:

``sampleRate``
^^^^^^^^^^^^^^

(Formerly ``sRateHz``)

Sampling rate (in Hz) of raw recording.

**Default** is 30000.

.. _shankMap:

``shankMap``
^^^^^^^^^^^^

(Formerly ``viShank_site``)

Shank ID of each site.

**Default** is empty.

.. _siteLoc:

``siteLoc``
^^^^^^^^^^^

(Formerly ``mrSiteXY``)

Site locations (in μm) (x values in the first column, y values in the second column).

**Default** is empty.

.. _siteMap:

``siteMap``
^^^^^^^^^^^

(Formerly ``viSite2Chan``)

Map of channel index to site ID (The mapping siteMap(i) = j corresponds to the statement 'site i is stored as channel j in the recording').

**Default** is empty.

.. _trialFile:

``trialFile``
^^^^^^^^^^^^^

(Formerly ``vcFile_trial``)

Path to file containing trial data (Can be .mat or .csv, must contain timestamps of trials in units of s).

**Default** is an empty string.

Advanced parameters
-------------------

.. _auxChan:

``auxChan``
^^^^^^^^^^^

(Formerly ``iChan_aux``)

Auxiliary channel index.

**Default** is empty.

.. _auxFile:

``auxFile``
^^^^^^^^^^^

(Formerly ``vcFile_aux``)

Path to file containing auxiliary channel.

**Default** is an empty string.

.. _auxLabel:

``auxLabel``
^^^^^^^^^^^^

(Formerly ``vcLabel_aux``)

Label for auxiliary channel data.

**Default** is 'Aux channel'.

.. _auxSampleRate:

``auxSampleRate``
^^^^^^^^^^^^^^^^^

(Formerly ``sRateHz_aux``)

Sample rate for auxiliary file.

**Default** is empty.

.. _auxScale:

``auxScale``
^^^^^^^^^^^^

(Formerly ``vrScale_aux``)

Scale factor for aux data.

**Default** is 1.

.. _batchMode:

``batchMode``
^^^^^^^^^^^^^

Suppress message boxes in favor of console messages.

**Default** is true.

.. _colorMap:

``colorMap``
^^^^^^^^^^^^

(Formerly ``mrColor_proj``)

RGB color map for background, primary selected, and secondary selected spikes (The first three values are the R values, the next three are the G values, and the last three are the B values.).

**Default** is [0.83203, 0, 0.9375, 0.85547, 0.50781, 0.46484, 0.91797, 0.76563, 0.085938].

.. _corrRange:

``corrRange``
^^^^^^^^^^^^^

(Formerly ``corrLim``)

Correlation score range to distinguish by color map.

**Default** is [0.9, 1].

.. _detectBipolar:

``detectBipolar``
^^^^^^^^^^^^^^^^^

(Formerly ``fDetectBipolar``)

Detect positive as well as negative peaks.

**Default** is false.

.. _dispFeature:

``dispFeature``
^^^^^^^^^^^^^^^

(Formerly ``vcFet_show``)

Feature to display in the feature projection plot.

One of the following:

- 'cov'
- 'pca'
- 'ppca'
- 'vpp'

**Default** is 'vpp'.

.. _dispFilter:

``dispFilter``
^^^^^^^^^^^^^^

(Formerly ``vcFilter_show``)

Filter to apply in traces plot.

One of the following:

- 'ndiff'
- 'sgdiff'
- 'bandpass'
- 'fir1'
- 'user'
- 'none'

**Default** is 'none'.

.. _driftMerge:

``driftMerge``
^^^^^^^^^^^^^^

(Formerly ``fDrift_merge``)

Compute multiple waveforms at three drift locations based on the spike position if true.

**Default** is true.

.. _evtManualThresh:

``evtManualThresh``
^^^^^^^^^^^^^^^^^^^

(Formerly ``spkThresh_uV``)

Manually-set spike detection threshold (in μV).

**Default** is empty.

.. _evtMergeRad:

``evtMergeRad``
^^^^^^^^^^^^^^^

(Formerly ``maxDist_site_um``)

Maximum distance (in μm) for merging spike waveforms.

**Default** is 50.

.. _evtWindowMergeFactor:

``evtWindowMergeFactor``
^^^^^^^^^^^^^^^^^^^^^^^^

(Formerly ``spkLim_factor_merge``)

Ratio of samples to take when computing correlation.

**Default** is 1.

.. _evtWindowRaw:

``evtWindowRaw``
^^^^^^^^^^^^^^^^

(Formerly ``spkLim_raw_ms``)

Time range (in ms) of raw spike waveforms, centered at the peak.

**Default** is [-0.5, 1.5].

.. _fftThresh:

``fftThresh``
^^^^^^^^^^^^^

(Formerly ``fft_thresh``)

Threshold (in MADs of power-frequency product) above which to remove frequency outliers.

**Default** is 0.

.. _figList:

``figList``
^^^^^^^^^^^

List of tags of figures to display in feature view.

One of the following:

- 'FigCorr'
- 'FigHist'
- 'FigISI'
- 'FigMap'
- 'FigPos'
- 'FigProj'
- 'FigRD'
- 'FigSim'
- 'FigTime'
- 'FigWav'

**Default** is ["FigCorr", "FigHist", "FigISI", "FigMap", "FigPos", "FigProj", "FigRD", "FigSim", "FigTime", "FigWav"].

.. _frFilterShape:

``frFilterShape``
^^^^^^^^^^^^^^^^^

(Formerly ``filter_shape_rate``)

Kernel shape for temporal averaging (Used in estimation of the firing rate of a given unit).

One of the following:

- 'triangle'
- 'rectangle'

**Default** is 'triangle'.

.. _frPeriod:

``frPeriod``
^^^^^^^^^^^^

(Formerly ``filter_sec_rate``)

Time period (in s) over which to determine firing rate (Used in estimation of the firing rate of a given unit).

**Default** is 2.

.. _frSampleRate:

``frSampleRate``
^^^^^^^^^^^^^^^^

(Formerly ``sRateHz_rate``)

Resampling rate (in Hz) for estimating the firing rate (Used in estimation of the firing rate of a given unit).

**Default** is 1000.

.. _freqLimNotch:

``freqLimNotch``
^^^^^^^^^^^^^^^^

Frequency ranges to exclude for notch filter.

**Default** is empty.

.. _freqLimStop:

``freqLimStop``
^^^^^^^^^^^^^^^

Frequency range to exclude for band-stop filter.

**Default** is empty.

.. _gainBoost:

``gainBoost``
^^^^^^^^^^^^^

(Formerly ``gain_boost``)

Scale factor to boost gain in raw recording (Used in filtering operation).

**Default** is 1.

.. _gpuLoadFactor:

``gpuLoadFactor``
^^^^^^^^^^^^^^^^^

GPU memory usage factor (Use 1/gpuLoadFactor amount of GPU memory).

**Default** is 5.

.. _groupShank:

``groupShank``
^^^^^^^^^^^^^^

(Formerly ``fGroup_shank``)

Group all sites on the same shank if true.

**Default** is true.

.. _gtFile:

``gtFile``
^^^^^^^^^^

(Formerly ``vcFile_gt``)

Path to file containing ground-truth data.

**Default** is an empty string.

.. _interpPC:

``interpPC``
^^^^^^^^^^^^

(Formerly ``fInterp_fet``)

Interpolate 1st principal vector to maximize projection of spikes if true.

**Default** is true.

.. _lfpSampleRate:

``lfpSampleRate``
^^^^^^^^^^^^^^^^^

(Formerly ``sRateHz_lfp``)

Sampling rate for LFP channels.

**Default** is 2500.

.. _loadTimeLimits:

``loadTimeLimits``
^^^^^^^^^^^^^^^^^^

(Formerly ``tlim_load``)

Time range (in s) of samples to load at once (All samples are loaded if empty).

**Default** is empty.

.. _maxAmp:

``maxAmp``
^^^^^^^^^^

Amplitude scale (in μV).

**Default** is 250.

.. _maxBytesLoad:

``maxBytesLoad``
^^^^^^^^^^^^^^^^

(Formerly ``MAX_BYTES_LOAD``)

Maximum number of bytes to load into memory.

**Default** is empty.

.. _maxClustersSite:

``maxClustersSite``
^^^^^^^^^^^^^^^^^^^

(Formerly ``maxCluPerSite``)

Maximum number of cluster centers computed per site (Used if RDDetrendMode is 'local').

**Default** is 20.

.. _maxSecLoad:

``maxSecLoad``
^^^^^^^^^^^^^^

(Formerly ``MAX_LOAD_SEC``)

Maximum sample duration (in s) to load into memory (Overrides maxBytesLoad if nonempty).

**Default** is empty.

.. _meanInterpFactor:

``meanInterpFactor``
^^^^^^^^^^^^^^^^^^^^

(Formerly ``nInterp_merge``)

Interpolation factor for mean unit waveforms (Set to 1 to disable).

**Default** is 1.

.. _minNeighborsDetect:

``minNeighborsDetect``
^^^^^^^^^^^^^^^^^^^^^^

(Formerly ``nneigh_min_detect``)

Minimum number of sample neighbors exceeding threshold for a sample to be considered a peak.

**Default** is 0.

.. _minSitesWeightFeatures:

``minSitesWeightFeatures``
^^^^^^^^^^^^^^^^^^^^^^^^^^

(Formerly ``min_sites_mask``)

Minimum number of sites to have if using weightFeatures (Ignored if weightFeatures is false).

**Default** is 5.

.. _nClustersShowAux:

``nClustersShowAux``
^^^^^^^^^^^^^^^^^^^^

(Formerly ``nClu_show_aux``)

Number of clusters to show in the aux vs. firing rate correlation.

**Default** is 10.

.. _nDiffOrder:

``nDiffOrder``
^^^^^^^^^^^^^^

(Formerly ``nDiff_filt``)

Order for differentiator filter (Used if and only if filterType is 'sgdiff' or 'ndiff').

**Default** is 2.

.. _nLoadsMaxPreview:

``nLoadsMaxPreview``
^^^^^^^^^^^^^^^^^^^^

(Formerly ``nLoads_max_preview``)

Number of time segments to load in preview.

**Default** is 30.

.. _nPassesMerge:

``nPassesMerge``
^^^^^^^^^^^^^^^^

(Formerly ``nRepeat_merge``)

Number of times to repeat automatic waveform-based merging.

**Default** is empty.

.. _nPeaksFeatures:

``nPeaksFeatures``
^^^^^^^^^^^^^^^^^^

(Formerly ``nFet_use``)

Number of potential peaks to use when computing features.

**Default** is 2.

.. _nSamplesPad:

``nSamplesPad``
^^^^^^^^^^^^^^^

(Formerly ``nPad_filt``)

Number of samples to overlap between chunks in large files.

**Default** is 100.

.. _nSecsLoadPreview:

``nSecsLoadPreview``
^^^^^^^^^^^^^^^^^^^^

(Formerly ``sec_per_load_preview``)

Number of seconds to load in preview.

**Default** is 1.

.. _nSegmentsTraces:

``nSegmentsTraces``
^^^^^^^^^^^^^^^^^^^

(Formerly ``nTime_traces``)

Number of time segments to display in traces view (A value of 1 shows one continuous time segment).

**Default** is 1.

.. _nSitesFigProj:

``nSitesFigProj``
^^^^^^^^^^^^^^^^^

Number of sites to show in feature projection view.

**Default** is 5.

.. _nSkip:

``nSkip``
^^^^^^^^^

(Formerly ``nSkip_show``)

Show every nSkip samples when plotting traces.

**Default** is 1.

.. _nSpikesFigISI:

``nSpikesFigISI``
^^^^^^^^^^^^^^^^^

Maximum number of spikes to show in ISI view.

**Default** is 200.

.. _nThreadsGPU:

``nThreadsGPU``
^^^^^^^^^^^^^^^

(Formerly ``nThreads``)

Number of GPU threads to use for clustering.

**Default** is 128.

.. _outlierThresh:

``outlierThresh``
^^^^^^^^^^^^^^^^^

(Formerly ``thresh_mad_clu``)

Threshold (in MADs) to remove outlier spikes for each cluster.

**Default** is 7.5.

.. _pcPair:

``pcPair``
^^^^^^^^^^

Pair of PCs to display.

**Default** is [1, 2].

.. _projTimeLimits:

``projTimeLimits``
^^^^^^^^^^^^^^^^^^

(Formerly ``tLimFigProj``)

Time range (in s) to display in feature projection view.

**Default** is empty.

.. _psthTimeBin:

``psthTimeBin``
^^^^^^^^^^^^^^^

(Formerly ``tbin_psth``)

Time bin (in s) for PSTH view.

**Default** is 0.01.

.. _psthXTick:

``psthXTick``
^^^^^^^^^^^^^

(Formerly ``xtick_psth``)

PSTH time tick mark spacing.

**Default** is 0.2.

.. _ramToGPUFactor:

``ramToGPUFactor``
^^^^^^^^^^^^^^^^^^

(Formerly ``nLoads_gpu``)

Ratio of RAM to GPU memory.

**Default** is 8.

.. _randomSeed:

``randomSeed``
^^^^^^^^^^^^^^

Seed for the random number generator.

**Default** is 0.

.. _showRaw:

``showRaw``
^^^^^^^^^^^

(Formerly ``fWav_raw_show``)

Show raw traces in waveform view if true.

**Default** is false.

.. _showSpikeCount:

``showSpikeCount``
^^^^^^^^^^^^^^^^^^

(Formerly ``fText``)

Show spike count per unit in waveform plot.

**Default** is true.

.. _siteCorrThresh:

``siteCorrThresh``
^^^^^^^^^^^^^^^^^^

(Formerly ``thresh_corr_bad_site``)

Threshold to reject bad sites based on maximum correlation with neighboring sites (Set to 0 to disable).

**Default** is 0.

.. _spikeThreshMax:

``spikeThreshMax``
^^^^^^^^^^^^^^^^^^

(Formerly ``spkThresh_max_uV``)

Maximum absolute amplitude (in μV) permitted for spikes.

**Default** is empty.

.. _tallSkinny:

``tallSkinny``
^^^^^^^^^^^^^^

(Formerly ``fTranspose_bin``)

Recording will be interpreted as nChannels x nSamples if true.

**Default** is true.

.. _threshFile:

``threshFile``
^^^^^^^^^^^^^^

(Formerly ``vcFile_thresh``)

Path to .mat file storing the spike detection threshold (Created by preview GUI).

**Default** is an empty string.

.. _umPerPix:

``umPerPix``
^^^^^^^^^^^^

(Formerly ``um_per_pix``)

Vertical site center-to-center spacing.

**Default** is 20.

.. _useElliptic:

``useElliptic``
^^^^^^^^^^^^^^^

(Formerly ``fEllip``)

Use elliptic (bandpass) filter if true (Uses Butterworth filter if false).

**Default** is true.

.. _useGPU:

``useGPU``
^^^^^^^^^^

(Formerly ``fGpu``)

Use GPU where appropriate.

**Default** is true.

.. _useGlobalDistCut:

``useGlobalDistCut``
^^^^^^^^^^^^^^^^^^^^

(Formerly ``fDc_global``)

Use a global distance cutoff for all sites if true.

**Default** is false.

.. _useParfor:

``useParfor``
^^^^^^^^^^^^^

(Formerly ``fParfor``)

Use parfor where appropriate.

**Default** is true.

.. _userFiltKernel:

``userFiltKernel``
^^^^^^^^^^^^^^^^^^

(Formerly ``vnFilter_user``)

User-specified filter kernel (Ignored unless filterType is 'user').

**Default** is empty.

.. _verbose:

``verbose``
^^^^^^^^^^^

(Formerly ``fVerbose``)

Be chatty when processing.

**Default** is true.

.. _weightFeatures:

``weightFeatures``
^^^^^^^^^^^^^^^^^^

(Formerly ``fSpatialMask_clu``)

Weight display features by distance from site if true.

**Default** is false.
