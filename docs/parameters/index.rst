.. _parameters:

JRCLUST parameters
==================

Common parameters
------------------

.. _CARMode:

CARMode
^^^^^^^

(Formerly ``vcCommonRef``)

The meaning of 'average' in 'common average reference'.

One of the following:

- 'mean': Computes the mean across sites, excluding sites in :ref:`ignoreSites`.
- 'median': Computes the median across sites, excluding sites in :ref:`ignoreSites`.
- 'none': Skips CAR.

**Default** is 'mean'.

.. _RDDetrendMode:

RDDetrendMode
^^^^^^^^^^^^^

(Formerly ``vcDetrend_postclu``)

Detrending mode to apply to rho-delta values in order to :ref:`determine cluster centers <assign-clusters>`.

One of the following:

- 'global': Performs the detrending for all sites together.
- 'local': Performs the detrending for each site separately.
- 'logz': Use z-scores instead of detrending.
- 'regress': Use a non-constant threshold to determine centers for each site.
- 'none': Skips detrending.

**Default** is 'global'.

.. _autoMergeBy:

autoMergeBy
^^^^^^^^^^^

(Formerly ``autoMergeCriterion``)

Metric to use for :ref:`automerging clusters <merge-post-hoc>` based on average waveform.

One of the following:

- 'pearson': Use the `Pearson correlation coefficient`_ to compute similarity.
- 'dist': Define similarity as :math:`1 - \frac{\|\text{tr}_i - \text{tr}_j\|}{\max(\|\text{tr}_i\|, \|\text{tr}_j\|)}`,
  where :math:`\text{tr}_k` is the mean waveform for cluster :math:`k`.

**Default** is 'pearson'.

.. _bitScaling:

bitScaling
^^^^^^^^^^

(Formerly ``uV_per_bit``)

ADC bit scaling factor (Conversion factor for ADC bit values to μV).

**Default** is 0.30518.

.. _blankPeriod:

blankPeriod
^^^^^^^^^^^

(Formerly ``blank_period_ms``)

Duration of blanking period (in ms) when the common mean exceeds :ref:`blankThresh`.

**Default** is 5.

.. _blankThresh:

blankThresh
^^^^^^^^^^^

(Formerly ``blank_thresh``)

Threshold (in MADs) above which to reject samples exceeding channel median after filtering.

**Default** is empty.

.. _clusterFeature:

clusterFeature
^^^^^^^^^^^^^^

(Formerly ``vcFet``)

The feature to extract from your spike waveforms in order to cluster them.

One of the following:

- 'cov': Covariance feature.
- 'energy': Standard deviation of spike waveform.
- 'gpca': :ref:`Global PCA <gpca>`.
- 'pca': Projection of spike waveforms onto :ref:`principal components <pca>`.
- 'vmin': :ref:`Minimum spike amplitude <vmin>`.
- 'vminmax': :ref:`Minimum and maximum spike amplitudes <vminmax>` (2 features/site).
- 'vpp': :ref:`Peak-to-peak spike amplitudes <feature-vpp>`.

**Default** is 'pca'.

.. _dataType:

dataType
^^^^^^^^

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

dispTimeLimits
^^^^^^^^^^^^^^

(Formerly ``tlim``)

Time range (in ms) to display.

**Default** is [0, 0.2].

.. _distCut:

distCut
^^^^^^^

(Formerly ``dc_percent``)

Percentile of pairwise distances between spikes on a site to use as a cutoff distance.

**Default** is 2.

.. _evtDetectRad:

evtDetectRad
^^^^^^^^^^^^

(Formerly ``maxDist_site_um``)

Maximum distance (in μm) to search over for :ref:`potential duplicates <merge-peaks>` (``r1`` in the figure below).
This distance is used to determine the number of sites to extract features if :ref:`nSiteDir` is empty.

**Default** is 50.

.. image:: /.static/evtDetectRad.png
   :scale: 25%

.. _evtGroupRad:

evtGroupRad
^^^^^^^^^^^^

(Formerly ``maxDist_site_spk_um``)

Maximum distance (in μm) for :ref:`extracting spike waveforms <extract-windows>`
(``r2`` in the figure below).

**Default** is 75.

.. image:: /.static/evtDetectRad.png
   :scale: 25%

.. _evtMergeRad:

evtMergeRad
^^^^^^^^^^^^

(Formerly ``maxDist_site_merge_um``)

Maximum distance (in μm) between sites to consider :ref:`merging a pair of units <merge-post-hoc>`.

**Default** is 35.

.. _evtWindow:

evtWindow
^^^^^^^^^

(Formerly ``spkLim_ms``)

Time range (in ms) of filtered spike waveforms, centered at the peak.

Must be an array with 2 elements, the first negative and the second positive.
For example, if ``evtWindow`` is set to [-0.5, 0.5], then 1/2 ms worth of samples
are extracted before and after the spiking event.

**Default** is [-0.25, 0.75].

.. _filtOrder:

filtOrder
^^^^^^^^^

Bandpass filter order.

**Default** is 3.

.. _filterType:

filterType
^^^^^^^^^^

(Formerly ``vcFilter``)

Type of filter to use on raw data.

One of the following:

- 'ndiff': Applies a differentiation filter, choosing a kernel depending on the order given in :ref:`nDiffOrder`.
- 'sgdiff': Applies a `Savitzky-Golay <https://en.wikipedia.org/wiki/Savitzky–Golay_filter>`_ filter depending on the order given in :ref:`nDiffOrder`.
- 'bandpass'
- 'fir1'
- 'user': Convolves your raw samples with a :ref:`kernel of your choosing <userFiltKernel>`.
- 'none': Skips filtering (not recommended).

**Default** is 'ndiff'.

.. _freqLimBP:

freqLimBP
^^^^^^^^^

(Formerly ``freqLim``)

Frequency cutoffs for bandpass filter.

**Default** is [300, 3000].

.. _headerOffset:

headerOffset
^^^^^^^^^^^^

(Formerly ``header_offset``)

Recording file header offset (in bytes).

JRCLUST will skip this many bytes at the beginning of your recording file.

**Default** is 0.

.. _ignoreChans:

ignoreChans
^^^^^^^^^^^

(Formerly ``viChanZero``)

Channel numbers to ignore manually.
Should be an array of integers between 1 and nChans.
These values will be merged into :ref:`ignoreSites`.

**Default** is empty.

.. _ignoreSites:

ignoreSites
^^^^^^^^^^^

(Formerly ``viSiteZero``)

Site IDs to ignore manually.
Should be an array of integers between 1 and nSites.

**Default** is empty.

.. _log10DeltaCut:

log10DeltaCut
^^^^^^^^^^^^^

(Formerly ``delta1_cut``)

:math:`\log_{10}` of delta cutoff (Spikes with delta values below this cutoff will not be considered as cluster centers).

**Default** is 0.6.

.. _log10RhoCut:

log10RhoCut
^^^^^^^^^^^

(Formerly ``rho_cut``)

:math:`\log_{10}` of rho cutoff (Spikes with rho values below this cutoff will not be considered as cluster centers).

**Default** is -2.5.

.. _maxUnitSim:

maxUnitSim
^^^^^^^^^^

(Formerly ``maxWavCor``)

Threshold for merging two units having similar spike waveforms (Units with a similiarity score above this value will be merged).

See :ref:`autoMergeBy` for how "similarity" is defined.

**Default** is 0.98.

.. _minClusterSize:

minClusterSize
^^^^^^^^^^^^^^

(Formerly ``min_count``)

Minimum number of spikes per cluster (Automatically set to the maximum of this value and twice the number of features).

**Default** is 30.

.. _nChans:

nChans
^^^^^^

Number of channels stored in recording file (Distinct from the number of AP sites).

**Default** is 384.

.. _nClusterIntervals:

nClusterIntervals
^^^^^^^^^^^^^^^^^

(Formerly ``nTime_clu``)

Number of intervals to divide the recording into around a spike.

When clustering, take the :math:`\frac{1}{\text{nClusterIntervals}}` fraction of all
spikes around a spiking event to compute distance.

For example, if ``nClusterIntervals`` = 1, all spikes will be used;
if ``nClusterIntervals`` = 2, JRCLUST will take the half of all spikes which are closest
in time to compute distances.
Increasing this value will take fewer and fewer spikes to compare at the risk of
oversplitting clusters (you might want to do this if you observe fast drift in your
recording).
However, automated merging based on the :ref:`waveform correlation <maxUnitSim>`
can merge most of the units initially split by drift.

**Default** is 4.

.. _nPCsPerSite:

nPCsPerSite
^^^^^^^^^^^

(Formerly ``nPcPerChan``)

Number of principal components to compute per site.

**Default** is 1.

.. _nSiteDir:

nSiteDir
^^^^^^^^

(Formerly ``maxSite``)

Number of neighboring sites to group in either direction.

The total number of sites per spike group (``nSitesEvt``) is 1 + 2\*``nSiteDir``.
In other words, a spike group includes the site on which the spike occurs, along with ``nSiteDir``
sites in the horizontal direction and ``nSiteDir`` in the vertical direction.

If empty, the number of sites per spike group is determined from :ref:`evtGroupRad`.

.. warning::
   This parameter may be deprecated in an upcoming release in favor of ``evtGroupRad``.

**Default** is empty.

.. _nSitesExcl:

nSitesExcl
^^^^^^^^^^

(Formerly ``nSites_ref``)

Number of sites to exclude from the spike waveform group for feature extraction.

**Default** is empty.

.. _nSpikesFigProj:

nSpikesFigProj
^^^^^^^^^^^^^^

(Formerly ``nShow_proj``)

Maximum number of spikes per cluster to display in the feature projection view.

**Default** is 500.

.. _nSpikesFigWav:

nSpikesFigWav
^^^^^^^^^^^^^

(Formerly ``nSpk_show``)

Maximum number of spikes per cluster to display generally.

**Default** is 30.

.. _outputDir:

outputDir
^^^^^^^^^^^^^

Directory in which to place output files (Will output to the same directory as this file if empty).

**Default** is an empty string.

.. _probePad:

probePad
^^^^^^^^^^^^

(Formerly ``vrSiteHW``)

Recording contact pad size (in μm) (Height x width).

**Default** is empty.

.. _psthTimeLimits:

psthTimeLimits
^^^^^^^^^^^^^^^^^^

(Formerly ``tlim_psth``)

Time range (in s) over which to display PSTH.

**Default** is empty.

.. _qqFactor:

qqFactor
^^^^^^^^^^^^

Spike detection threshold.

Multiplier of the :ref:`estimate <compute-threshold>` :math:`\sigma_{\text{noise}}^{(i)}`
of standard deviation of noise distribution on each site to compute the threshold for that site.
In other words,

.. math::

    \text{Thr}_i := \text{qqFactor} \cdot \sigma_{\text{noise}}^{(i)}

is the spike detection threshold for site :math:`i`.

**Default** is 5.

.. _rawRecordings:

rawRecordings
^^^^^^^^^^^^^^^^^

Path or paths to raw recordings to sort.

**Default** is empty.

.. _recordingFormat:

recordingFormat
^^^^^^^^^^^^^^^

Format of raw recordings.

One of the following:

- 'SpikeGLX': A flat [binary file])(https://github.com/billkarsh/SpikeGLX/blob/master/Markdown/UserManual.md#output-file-format) with no header, stored in channels x samples order.
- 'Intan': Traditional [Intan file format](http://www.intantech.com/files/Intan_RHD2000_data_file_formats.pdf).

**Default** is 'SpikeGLX'.

.. _refracInt:

refracInt
^^^^^^^^^

(Formerly ``spkRefrac_ms``)

Spike refractory period (in ms).

**Default** is 0.25.

.. _sampleRate:

sampleRate
^^^^^^^^^^

(Formerly ``sRateHz``)

Sampling rate (in Hz) of raw recording.

**Default** is 30000.

.. _shankMap:

shankMap
^^^^^^^^^^^^

(Formerly ``viShank_site``)

Shank ID of each site.

**Default** is empty.

.. _siteLoc:

siteLoc
^^^^^^^^^^^

(Formerly ``mrSiteXY``)

Site locations (in μm) (x values in the first column, y values in the second column).

**Default** is empty.

.. _siteMap:

siteMap
^^^^^^^^^^^

(Formerly ``viSite2Chan``)

Map of channel index to site ID (The mapping siteMap(i) = j corresponds to the statement 'site i is stored as channel j in the recording').

**Default** is empty.

.. _trialFile:

trialFile
^^^^^^^^^^^^^

(Formerly ``vcFile_trial``)

Path to file containing trial data (Can be .mat or .csv, must contain timestamps of trials in units of s).

**Default** is an empty string.

Advanced parameters
-------------------

.. _auxChan:

auxChan
^^^^^^^^^^^

(Formerly ``iChan_aux``)

Auxiliary channel index.

**Default** is empty.

.. _auxFile:

auxFile
^^^^^^^^^^^

(Formerly ``vcFile_aux``)

Path to file containing auxiliary channel.

**Default** is an empty string.

.. _auxLabel:

auxLabel
^^^^^^^^^^^^

(Formerly ``vcLabel_aux``)

Label for auxiliary channel data.

**Default** is 'Aux channel'.

.. _auxSampleRate:

auxSampleRate
^^^^^^^^^^^^^^^^^

(Formerly ``sRateHz_aux``)

Sample rate for auxiliary file.

**Default** is empty.

.. _auxScale:

auxScale
^^^^^^^^^^^^

(Formerly ``vrScale_aux``)

Scale factor for aux data.

**Default** is 1.

.. _batchMode:

batchMode
^^^^^^^^^^^^^

Suppress message boxes in favor of console messages.

**Default** is true.

.. _colorMap:

colorMap
^^^^^^^^^^^^

(Formerly ``mrColor_proj``)

RGB color map for background, primary selected, and secondary selected spikes (The first three values are the R values, the next three are the G values, and the last three are the B values.).

**Default** is [0.83203, 0, 0.9375, 0.85547, 0.50781, 0.46484, 0.91797, 0.76563, 0.085938].

.. _corrRange:

corrRange
^^^^^^^^^^^^^

(Formerly ``corrLim``)

Correlation score range to distinguish by color map.

**Default** is [0.9, 1].

.. _detectBipolar:

detectBipolar
^^^^^^^^^^^^^^^^^

(Formerly ``fDetectBipolar``)

Detect positive as well as negative peaks.

**Default** is false.

.. _dispFeature:

dispFeature
^^^^^^^^^^^^^^^

(Formerly ``vcFet_show``)

Feature to display in the feature projection plot.

One of the following:

- 'cov'
- 'pca'
- 'ppca'
- 'vpp'
- 'template' (only supported for template-based clusterings, e.g., Kilosort)

**Default** is 'vpp'.

.. _dispFilter:

dispFilter
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

driftMerge
^^^^^^^^^^^^^^

(Formerly ``fDrift_merge``)

Compute multiple waveforms at three drift locations based on the spike position if true.

**Default** is true.

.. _evtManualThresh:

evtManualThresh
^^^^^^^^^^^^^^^

(Formerly ``spkThresh_uV``)

Manually-set spike detection threshold (in μV).

**Default** is empty.

.. _evtWindowMergeFactor:

evtWindowMergeFactor
^^^^^^^^^^^^^^^^^^^^^^^^

(Formerly ``spkLim_factor_merge``)

Ratio of samples to take when computing correlation.

**Default** is 1.

.. _evtWindowRaw:

evtWindowRaw
^^^^^^^^^^^^^^^^

(Formerly ``spkLim_raw_ms``)

Time range (in ms) of raw spike waveforms, centered at the peak.

Must be an array with 2 elements, the first negative and the second positive.
For example, if ``evtWindowRaw`` is set to [-1, 1], then 1 ms worth of samples
are extracted before and after the spiking event.

**Default** is [-0.5, 1.5].

.. _extractAfterDetect:

extractAfterDetect
^^^^^^^^^^^^^^^^^^

Extract features only after detecting all spikes across all files if true.
Otherwise, features will be computed on a per-chunk, per-file basis, which may not be what you want.

This is effectively set to true if you specify :ref:`clusterFeature` = 'gpca'.

**Default** is false.

.. _fftThresh:

fftThresh
^^^^^^^^^^^^^

(Formerly ``fft_thresh``)

Threshold (in MADs of power-frequency product) above which to remove frequency outliers
when :ref:`denoising <denoising>`.
Frequencies with power-frequency product above this threshold will be zeroed out as noise.

Setting to 0 disables this notch filtering.
If you choose to enable, the recommended value is 10.

**Default** is 0.

.. _figList:

figList
^^^^^^^

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

frFilterShape
^^^^^^^^^^^^^^^^^

(Formerly ``filter_shape_rate``)

Kernel shape for temporal averaging (Used in estimation of the firing rate of a given unit).

One of the following:

- 'triangle'
- 'rectangle'

**Default** is 'triangle'.

.. _frPeriod:

frPeriod
^^^^^^^^^^^^

(Formerly ``filter_sec_rate``)

Time period (in s) over which to determine firing rate (Used in estimation of the firing rate of a given unit).

**Default** is 2.

.. _frSampleRate:

frSampleRate
^^^^^^^^^^^^^^^^

(Formerly ``sRateHz_rate``)

Resampling rate (in Hz) for estimating the firing rate (Used in estimation of the firing rate of a given unit).

**Default** is 1000.

.. _freqLimNotch:

freqLimNotch
^^^^^^^^^^^^^^^^

Frequency ranges to exclude for notch filter.

**Default** is empty.

.. _freqLimStop:

freqLimStop
^^^^^^^^^^^^^^^

Frequency range to exclude for band-stop filter.

**Default** is empty.

.. _gainBoost:

gainBoost
^^^^^^^^^^^^^

(Formerly ``gain_boost``)

Scale factor to boost gain in raw recording (Used in filtering operation).

**Default** is 1.

.. _gpuLoadFactor:

gpuLoadFactor
^^^^^^^^^^^^^^^^^

GPU memory usage factor (Use 1/gpuLoadFactor amount of GPU memory).

**Default** is 5.

.. _groupShank:

groupShank
^^^^^^^^^^^^^^

(Formerly ``fGroup_shank``)

Group all sites on the same shank if true.

**Default** is true.

.. _gtFile:

gtFile
^^^^^^^^^^

(Formerly ``vcFile_gt``)

Path to file containing ground-truth data.

**Default** is an empty string.

.. _interpPC:

interpPC
^^^^^^^^

(Formerly ``fInterp_fet``)

Interpolate 1st principal vector to maximize projection of spikes if true.

**Default** is true.

.. _lfpSampleRate:

lfpSampleRate
^^^^^^^^^^^^^^^^^

(Formerly ``sRateHz_lfp``)

Sampling rate for LFP channels.

**Default** is 2500.

.. _loadTimeLimits:

loadTimeLimits
^^^^^^^^^^^^^^^^^^

(Formerly ``tlim_load``)

Time range (in s) of samples to load at once (All samples are loaded if empty).

**Default** is empty.

.. _maxAmp:

maxAmp
^^^^^^^^^^

Amplitude scale (in μV).

**Default** is 250.

.. _maxBytesLoad:

maxBytesLoad
^^^^^^^^^^^^^^^^

(Formerly ``MAX_BYTES_LOAD``)

Maximum number of bytes to load into memory.

**Default** is empty.

.. _maxClustersSite:

maxClustersSite
^^^^^^^^^^^^^^^^^^^

(Formerly ``maxCluPerSite``)

Maximum number of cluster centers computed per site (Used if :ref:`RDDetrendMode` is 'local').

**Default** is 20.

.. _maxSecLoad:

maxSecLoad
^^^^^^^^^^^^^^

(Formerly ``MAX_LOAD_SEC``)

Maximum sample duration (in s) to load into memory (Overrides maxBytesLoad if nonempty).

**Default** is empty.

.. _meanInterpFactor:

meanInterpFactor
^^^^^^^^^^^^^^^^^^^^

(Formerly ``nInterp_merge``)

Interpolation factor for mean unit waveforms (Set to 1 to disable).

**Default** is 1.

.. _minNeighborsDetect:

minNeighborsDetect
^^^^^^^^^^^^^^^^^^^^^^

(Formerly ``nneigh_min_detect``)

Minimum number of sample neighbors exceeding threshold for a sample to be considered a peak.

For example, consider a potential peak occurring at sample :math:`t_i`.
If ``minNeighborsDetect`` is set to 1, then **either** sample :math:`t_{i-1}` or :math:`t_{i+1}`
must also exceed the detection threshold.
If ``minNeighborsDetect`` is set to 2, then **both** sample :math:`t_{i-1}` and :math:`t_{i+1}`
must also exceed the detection threshold.
If ``minNeighborsDetect`` is set to 0, then samples :math:`t_{i-1}` and :math:`t_{i+1}`
are not considered.

Must be one of 0, 1, or 2.

**Default** is 0.

.. _minSitesWeightFeatures:

minSitesWeightFeatures
^^^^^^^^^^^^^^^^^^^^^^^^^^

(Formerly ``min_sites_mask``)

Minimum number of sites to have if using weightFeatures (Ignored if weightFeatures is false).

**Default** is 5.

.. _nClustersShowAux:

nClustersShowAux
^^^^^^^^^^^^^^^^^^^^

(Formerly ``nClu_show_aux``)

Number of clusters to show in the aux vs. firing rate correlation.

**Default** is 10.

.. _nDiffOrder:

nDiffOrder
^^^^^^^^^^^^^^

(Formerly ``nDiff_filt``)

Order for differentiator filter (Used if and only if :ref:`filterType` is 'sgdiff' or 'ndiff').

**Default** is 2.

.. _nLoadsMaxPreview:

nLoadsMaxPreview
^^^^^^^^^^^^^^^^^^^^

(Formerly ``nLoads_max_preview``)

Number of time segments to load in preview.

**Default** is 30.

.. _nPassesMerge:

nPassesMerge
^^^^^^^^^^^^^^^^

(Formerly ``nRepeat_merge``)

Number of times to repeat automatic waveform-based merging.

**Default** is 10.

.. _nPeaksFeatures:

nPeaksFeatures
^^^^^^^^^^^^^^^^^^

(Formerly ``nFet_use``)

Number of potential peaks to use when computing features.

**Default** is 2.

.. _nSamplesPad:

nSamplesPad
^^^^^^^^^^^^^^^

(Formerly ``nPad_filt``)

Number of samples to overlap between chunks in large files.

**Default** is 100.

.. _nSecsLoadPreview:

nSecsLoadPreview
^^^^^^^^^^^^^^^^^^^^

(Formerly ``sec_per_load_preview``)

Number of seconds to load in preview.

**Default** is 1.

.. _nSegmentsTraces:

nSegmentsTraces
^^^^^^^^^^^^^^^^^^^

(Formerly ``nTime_traces``)

Number of time segments to display in traces view (A value of 1 shows one continuous time segment).

**Default** is 1.

.. _nSitesFigProj:

nSitesFigProj
^^^^^^^^^^^^^^^^^

Number of sites to show in feature projection view.

**Default** is 5.

.. _nSkip:

nSkip
^^^^^^^^^

(Formerly ``nSkip_show``)

Show every nSkip samples when plotting traces.

**Default** is 1.

.. _nSpikesFigISI:

nSpikesFigISI
^^^^^^^^^^^^^^^^^

Maximum number of spikes to show in ISI view.

**Default** is 200.

.. _nThreadsGPU:

nThreadsGPU
^^^^^^^^^^^^^^^

(Formerly ``nThreads``)

Number of GPU threads to use for clustering.

**Default** is 128.

.. _outlierThresh:

outlierThresh
^^^^^^^^^^^^^^^^^

(Formerly ``thresh_mad_clu``)

Threshold (in MADs) to remove outlier spikes for each cluster.

**Default** is 7.5.

.. _pcPair:

pcPair
^^^^^^^^^^

Pair of PCs to display.

**Default** is [1, 2].

.. _projTimeLimits:

projTimeLimits
^^^^^^^^^^^^^^^^^^

(Formerly ``tLimFigProj``)

Time range (in s) to display in feature projection view.

**Default** is empty.

.. _psthTimeBin:

psthTimeBin
^^^^^^^^^^^^^^^

(Formerly ``tbin_psth``)

Time bin (in s) for PSTH view.

**Default** is 0.01.

.. _psthXTick:

psthXTick
^^^^^^^^^^^^^

(Formerly ``xtick_psth``)

PSTH time tick mark spacing.

**Default** is 0.2.

.. _ramToGPUFactor:

ramToGPUFactor
^^^^^^^^^^^^^^^^^^

(Formerly ``nLoads_gpu``)

Ratio of RAM to GPU memory.

**Default** is 8.

.. _randomSeed:

randomSeed
^^^^^^^^^^

Seed for the random number generator.

**Default** is 0.

.. _realignTraces:

realignTraces
^^^^^^^^^^^^^

Realign spike traces after subtracting local CAR (Realign if 1, perform subpixel interpolation if 2).

**Default** is 0.

.. _showRaw:

showRaw
^^^^^^^

(Formerly ``fWav_raw_show``)

Show raw traces in waveform view if true.

**Default** is false.

.. _showSpikeCount:

showSpikeCount
^^^^^^^^^^^^^^^^^^

(Formerly ``fText``)

Show spike count per unit in waveform plot.

**Default** is true.

.. _siteCorrThresh:

siteCorrThresh
^^^^^^^^^^^^^^^^^^

(Formerly ``thresh_corr_bad_site``)

Threshold to reject bad sites based on maximum correlation with neighboring sites (Set to 0 to disable).

**Default** is 0.

.. _spikeThreshMax:

spikeThreshMax
^^^^^^^^^^^^^^^^^^

(Formerly ``spkThresh_max_uV``)

Maximum absolute amplitude (in μV) permitted for spikes.

**Default** is empty.

.. _tallSkinny:

tallSkinny
^^^^^^^^^^^^^^

(Formerly ``fTranspose_bin``)

Recording will be interpreted as nChannels x nSamples if true.

**Default** is true.

.. _threshFile:

threshFile
^^^^^^^^^^^^^^

(Formerly ``vcFile_thresh``)

Path to .mat file storing the spike detection threshold (Created by preview GUI).

**Default** is an empty string.

.. _umPerPix:

umPerPix
^^^^^^^^^^^^

(Formerly ``um_per_pix``)

Vertical site center-to-center spacing.

**Default** is 20.

.. _useElliptic:

useElliptic
^^^^^^^^^^^^^^^

(Formerly ``fEllip``)

Use elliptic (bandpass) filter if true (Uses Butterworth filter if false).

**Default** is true.

.. _useGPU:

useGPU
^^^^^^^^^^

(Formerly ``fGpu``)

Use GPU where appropriate.

**Default** is true.

.. _useGlobalDistCut:

useGlobalDistCut
^^^^^^^^^^^^^^^^^^^^

(Formerly ``fDc_global``)

Use a global distance cutoff for all sites if true.

**Default** is false.

.. _useParfor:

useParfor
^^^^^^^^^^^^^

(Formerly ``fParfor``)

Use parfor where appropriate.

**Default** is true.

.. _userFiltKernel:

userFiltKernel
^^^^^^^^^^^^^^^^^^

(Formerly ``vnFilter_user``)

User-specified filter kernel (Ignored unless :ref:`filterType` is 'user').

Your filtered samples will be the output of a convolution of your raw samples
with this kernel.
You must specify this if and only if your :ref:`filterType` is ``'user'``.

**Default** is empty.

.. _verbose:

verbose
^^^^^^^^^^^

(Formerly ``fVerbose``)

Be chatty when processing.

**Default** is true.

.. _weightFeatures:

weightFeatures
^^^^^^^^^^^^^^^^^^

(Formerly ``fSpatialMask_clu``)

Weight display features by distance from site if true.

**Default** is false.

.. _`Pearson correlation coefficient`: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
.. _`MAD`: https://en.wikipedia.org/wiki/Median_absolute_deviation
