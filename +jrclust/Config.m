classdef Config < dynamicprops
    %CONFIG JRCLUST session configuration
    % replacement for P struct

    %% OBJECT-LEVEL PROPERTIES
    properties (Hidden, SetAccess=private, SetObservable)
        configFile;
        errMsg;
        isLoaded;
        isError;
        oldPcount;
        paramSet;
        probeFile;
        oldParamSet;
        tempParams;
    end

    %% OLD-STYLE PARAMS, publicly settable (will be deprecated after a grace period)
    properties (SetObservable, Dependent, Hidden, Transient)
        MAX_BYTES_LOAD; % => maxBytesLoad
        MAX_LOAD_SEC; % => maxSecLoad
        autoMergeCriterion; % => autoMergeBy
        blank_period_ms; % => blankPeriod
        blank_thresh; % => blankThresh
        corrLim; % => corrRange
        csFile_merge; % => multiRaw
        dc_percent; % => distCut
        delta1_cut; % => log10DeltaCut
        fDc_global; % => useGlobalDistCut
        fDetectBipolar; % => detectBipolar
        fDrift_merge; % => driftMerge
        fEllip; % => useElliptic
        fGpu; % => useGPU
        fGroup_shank; % => groupShank
        fInterp_fet; % => interpPC
        fParfor; % => useParfor
        fSpatialMask_clu; % => weightFeatures
        fText; % => showSpikeCount
        fTranspose_bin; % => tallSkinny
        fVerbose; % => verbose
        fWav_raw_show; % => showRaw
        fft_thresh; % => fftThreshMad
        filter_sec_rate; % => frPeriod
        filter_shape_rate; % => frFilterShape
        freqLim; % => freqLimBP
        gain_boost; % => gainBoost
        header_offset; % => headerOffset
        iChan_aux; % => auxChan
        maxCluPerSite; % => maxClustersSite
        maxDist_site_spk_um; % => evtDetectRad
        maxDist_site_um; % => evtMergeRad
        maxSite; % => nSiteDir
        maxWavCor; % => maxUnitSim
        min_count; % => minClusterSize
        min_sites_mask; % => minSitesWeightFeatures
        mrColor_proj; % => colorMap
        mrSiteXY; % => siteLoc
        nClu_show_aux; % => nClustersShowAux
        nDiff_filt; % => nDiffOrder
        nFet_use; % => nPeaksFeatures
        nInterp_merge; % => meanInterpFactor
        nLoads_gpu; % => ramToGPUFactor
        nLoads_max_preview; % => nLoadsMaxPreview
        nPad_filt; % => nSamplesPad
        nPcPerChan; % => nPCsPerSite
        nRepeat_merge; % => nPassesMerge
        nShow_proj; % => nSpikesFigProj
        nSites_ref; % => nSitesExcl
        nSkip_lfp; % => lfpDsFactor
        nSkip_show; % => nSkip
        nSpk_show; % => nSpikesFigWav
        nThreads; % => nThreadsGPU
        nTime_clu; % => nClusterIntervals
        nTime_traces; % => nSegmentsTraces
        nneigh_min_detect; % => minNeighborsDetect
        probe_file; % => probeFile
        rho_cut; % => log10RhoCut
        sRateHz; % => sampleRate
        sRateHz_aux; % => auxSampleRate
        sRateHz_lfp; % => lfpSampleRate
        sRateHz_rate; % => frSampleRate
        sec_per_load_preview; % => nSecsLoadPreview
        spkLim; % => evtWindowSamp
        spkLim_factor_merge; % => evtWindowMergeFactor
        spkLim_ms; % => evtWindow
        spkLim_raw; % => evtWindowRawSamp
        spkLim_raw_ms; % => evtWindowRaw
        spkRefrac; % => refracIntSamp
        spkRefrac_ms; % => refracInt
        spkThresh; % => evtManualThreshSamp
        spkThresh_max_uV; % => spikeThreshMax
        spkThresh_uV; % => evtManualThresh
        tLimFigProj; % => projTimeLimits
        tbin_psth; % => psthTimeBin
        thresh_corr_bad_site; % => siteCorrThresh
        thresh_mad_clu; % => outlierThresh
        time_feature_factor; % => timeFeatureFactor
        tlim; % => dispTimeLimits
        tlim_load; % => loadTimeLimits
        tlim_psth; % => psthTimeLimits
        uV_per_bit; % => bitScaling
        um_per_pix; % => umPerPix
        vcCommonRef; % => CARMode
        vcDataType; % => dataType
        vcDetrend_postclu; % => RDDetrendMode
        vcFet; % => clusterFeature
        vcFet_show; % => dispFeature
        vcFile; % => singleRaw
        vcFile_aux; % => auxFile
        vcFile_gt; % => gtFile
        vcFile_prm; % => configFile
        vcFile_thresh; % => threshFile
        vcFile_trial; % => trialFile
        vcFilter; % => filterType
        vcFilter_show; % => dispFilter
        vcLabel_aux; % => auxLabel
        viShank_site; % => shankMap
        viSite2Chan; % => siteMap
        viSiteZero; % => ignoreSites
        vnFilter_user; % => userFiltKernel
        vrScale_aux; % => auxScale
        vrSiteHW; % => probePad
        xtick_psth; % => psthXTick
    end

    %% OLD-STLYE PARAMS, not publicly settable
    properties (Dependent, SetAccess=private, Hidden)
        miSites;                    % => siteNeighbors
    end

    %% NEW-STYLE PARAMS, publicly settable
    properties (SetObservable)
        % USAGE PARAMETERS
        batchMode = true; % Suppress message boxes in favor of console messages
        outputDir = ''; % Directory in which to place output files
        verbose = true; % Be chatty when processing

        % EXECUTION PARAMETERS
        gpuLoadFactor = 5; % GPU memory usage factor
        maxBytesLoad = []; % Maximum number of bytes to load into memory
        maxSecLoad = []; % Maximum sample duration (in s) to load into memory
        nThreadsGPU = 128; % Number of GPU threads to use for clustering
        ramToGPUFactor = 8; % Ratio of RAM to GPU memory
        randomSeed = 0; % Seed for the random number generator
        useGPU = true; % Use GPU where appropriate
        useParfor = true; % Use parfor where appropriate

        % PROBE PARAMETERS
        probePad = []; % Recording contact pad size (in ?m)
        shankMap = []; % Shank ID of each site
        siteLoc = []; % Site locations (in ?m)
        siteMap = []; % Map of channel index to site ID

        % RECORDING FILE PARAMETERS
        bitScaling = 0.30518; % ADC bit scaling factor
        dataType = 'int16'; % Format of raw recordings
        headerOffset = 0; % Recording file header offset (in bytes)
        nChans = 384; % Number of channels stored in recording file
        rawRecordings = {}; % Path or paths to raw recordings to sort
        sampleRate = 30000; % Sampling rate (in Hz) of raw recording
        tallSkinny = true; % Recording will be interpreted as nChannels x nSamples if true

        % PREPROCESSING PARAMETERS
        ignoreSites = []; % Sites to ignore manually
        loadTimeLimits = []; % Time range (in s) of samples to load at once
        nDiffOrder = 2; % Order for differentiator filter
        nSamplesPad = 100; % Number of samples to overlap between chunks in large files
        useElliptic = true; % Use elliptic (bandpass) filter if true
        userFiltKernel = []; % User-specified filter kernel

        % SPIKE DETECTION PARAMETERS
        CARMode = 'mean'; % The meaning of 'average' in 'common average reference'
        blankPeriod = 5; % Duration of blanking period (in ms) when the common mean exceeds blankThresh
        blankThresh = []; % Threshold (in MADs) above which to reject samples exceeding channel median after filtering
        detectBipolar = false; % Detect positive as well as negative peaks
        evtDetectRad = 75; % Maximum distance (in ?m) for extracting spike waveforms
        evtManualThresh = []; % Manually-set spike detection threshold (in ?V)
        evtMergeRad = 50; % Maximum distance (in ?m) for merging spike waveforms
        evtWindow = [-0.25, 0.75]; % Time range (in ms) of filtered spike waveforms, centered at the peak
        evtWindowRaw = [-0.5, 1.5]; % Time range (in ms) of raw spike waveforms, centered at the peak
        fftThresh = 0; % Threshold (in MADs of power-frequency product) above which to remove frequency outliers
        filtOrder = 3; % Bandpass filter order
        filterType = 'ndiff'; % Type of filter to use on raw data
        freqLimBP = [300, 3000]; % Frequency cutoffs for bandpass filter
        freqLimNotch = []; % Frequency ranges to exclude for notch filter
        freqLimStop = []; % Frequency range to exclude for band-stop filter
        gainBoost = 1; % Scale factor to boost gain in raw recording
        groupShank = true; % Group all sites on the same shank if true
        minNeighborsDetect = 0; % Minimum number of sample neighbors exceeding threshold for a sample to be considered a peak
        nSiteDir = []; % Number of neighboring sites to group in either direction
        nSitesExcl = []; % Number of sites to exclude from the spike waveform group
        qqFactor = 5; % Spike detection threshold
        refracInt = 0.25; % Spike refractory period (in ms)
        spikeThreshMax = []; % Maximum absolute amplitude (in ?V) permitted for spikes
        threshFile = ''; % Path to .mat file storing the spike detection threshold

        % FEATURE EXTRACTION PARAMETERS
        clusterFeature = 'pca'; % The feature to extract from your spike waveforms in order to cluster them
        interpPC = true; % Interpolate 1st principal vector to maximize projection of spikes if true
        nPCsPerSite = 1; % Number of principal components to compute per site
        nPeaksFeatures = 2; % Number of potential peaks to use when computing features

        % CLUSTERING PARAMETERS
        RDDetrendMode = 'global'; % Detrending mode to apply to rho-delta values in order to determine cluster centers
        autoMergeBy = 'pearson'; % Metric to use for automerging clusters
        distCut = 2; % Percentile of pairwise distances between spikes on a site to use as a cutoff distance
        driftMerge = true; % Compute multiple waveforms at three drift locations based on the spike position if true
        evtWindowMergeFactor = 1; % Ratio of samples to take when computing correlation
        log10DeltaCut = 0.6; % Log10 of delta cutoff
        log10RhoCut = -2.5; % Log10 of rho cutoff
        maxClustersSite = 20; % Maximum number of cluster centers computed per site
        maxUnitSim = 0.98; % Threshold for merging two units having similar spike waveforms
        meanInterpFactor = 1; % Interpolation factor for mean unit waveforms
        minClusterSize = 30; % Minimum number of spikes per cluster
        minSitesWeightFeatures = 5; % Minimum number of sites to have if using weightFeatures
        nClusterIntervals = 4; % Number of intervals to divide the recording into around a spike
        nPassesMerge = []; % Number of times to repeat automatic waveform-based merging
        outlierThresh = 7.5; % Threshold (in MADs) to remove outlier spikes for each cluster
        useGlobalDistCut = false; % Use a global distance cutoff for all sites if true
        weightFeatures = false; % Weight display features by distance from site if true

        % CURATION PARAMETERS
        figList = ["FigCorr", "FigHist", "FigISI", "FigMap", "FigPos", "FigProj", "FigRD", "FigSim", "FigTime", "FigWav"]; % List of tags of figures to display in feature view
        frFilterShape = 'triangle'; % Kernel shape for temporal averaging
        frPeriod = 2; % Time period (in s) over which to determine firing rate
        frSampleRate = 1000; % Resampling rate (in Hz) for estimating the firing rate

        % DISPLAY PARAMETERS
        colorMap = [0.83203, 0, 0.9375, 0.85547, 0.50781, 0.46484, 0.91797, 0.76563, 0.085938]; % RGB color map for background, primary selected, and secondary selected spikes
        corrRange = [0.9, 1]; % Correlation score range to distinguish by color map
        dispFeature = 'vpp'; % Feature to display in the feature projection plot
        dispFilter = 'none'; % Filter to apply in traces plot
        dispTimeLimits = [0, 0.2]; % Time range (in ms) to display
        maxAmp = 250; % Amplitude scale (in ?V)
        nSitesFigProj = 5; % Number of sites to show in feature projection view
        nSpikesFigISI = 200; % Maximum number of spikes to show in ISI view
        nSpikesFigProj = 500; % Maximum number of spikes per cluster to display in the feature projection view
        nSpikesFigWav = 30; % Maximum number of spikes per cluster to display generally
        pcPair = [1, 2]; % Pair of PCs to display
        projTimeLimits = []; % Time range (in s) to display in feature projection view
        showRaw = false; % Show raw traces in waveform view if true
        showSpikeCount = true; % Show spike count per unit in waveform plot
        umPerPix = 20; % Vertical site center-to-center spacing

        % TRIAL PARAMETERS
        psthTimeBin = 0.01; % Time bin (in s) for PSTH view
        psthTimeLimits = []; % Time range (in s) over which to display PSTH
        psthXTick = 0.2; % PSTH time tick mark spacing
        trialFile = ''; % Path to file containing trial data

        % VALIDATION PARAMETERS
        gtFile = ''; % Path to file containing ground-truth data

        % PREVIEW PARAMETERS
        nLoadsMaxPreview = 30; % Number of time segments to load in preview
        nSecsLoadPreview = 1; % Number of seconds to load in preview
        siteCorrThresh = 0; % Threshold to reject bad sites based on maximum correlation with neighboring sites

        % TRACES PARAMETERS
        nSegmentsTraces = 1; % Number of time segments to display in traces view
        nSkip = 1; % Show every nSkip samples when plotting traces

        % LFP PARAMETERS
        lfpSampleRate = 2500; % Sampling rate for LFP channels

        % AUX CHANNEL PARAMETERS
        auxChan = []; % Auxiliary channel index
        auxFile = ''; % Path to file containing auxiliary channel
        auxLabel = 'Aux channel'; % Label for auxiliary channel data
        auxSampleRate = []; % Sample rate for auxiliary file
        auxScale = 1; % Scale factor for aux data
        nClustersShowAux = 10; % Number of clusters to show in the aux vs. firing rate correlation
    end

    %% NEW-STYLE PARAMS, not publicly settable
    properties (SetAccess=private)
        siteNeighbors;              % indices of neighbors for each site
    end

    %% COMPUTED PARAMS
    properties (SetObservable, Dependent)
        bytesPerSample;             % byte count for each discrete sample
        evtManualThreshSamp;        % evtManualThresh / bitScaling
        evtWindowRawSamp;           % interval around event to extract raw spike waveforms, in samples
        evtWindowSamp;              % interval around event to extract filtered spike waveforms, in samples
        nSites;                     % numel(siteMap)
        nSitesEvt;                  % 2*nSiteDir + 1 - nSitesExcl
        refracIntSamp;              % spike refractory interval, in samples
        sessionName;                % name of prm file, without path or extensions
    end

    %% LIFECYCLE
    methods
        function obj = Config(filename)
            %CONFIG Construct an instance of this class
            if nargin == 0 || isempty(filename)
                success = obj.makeCfg();
                if ~success
                    return;
                end

                filename = obj.configFile;
            end

            % read in default parameter set
            fid = fopen(fullfile(jrclust.utils.basedir(), 'params.json'), 'r');
            obj.paramSet = jsondecode(fread(fid, '*char'));
            fclose(fid);

            % read in mapping to old (v3) parameter set
            fid = fopen(fullfile(jrclust.utils.basedir(), 'old2new.json'), 'r');
            obj.oldParamSet = jsondecode(fread(fid, '*char'));
            fclose(fid);

            obj.isLoaded = 0;
            obj.isError = 0;

            obj.tempParams = containers.Map();

            if ~isfile(filename)
                obj.error(sprintf('Cannot load config: file %s not found', filename), 'Missing config file');
                return;
            elseif isempty(which(filename))
                obj.configFile = filename;
            else
                obj.configFile = which(filename);
            end

            obj.loadParams();
            if ~obj.isError
                obj.validateParams();
            end
            obj.isLoaded = ~obj.isError;
        end
    end

    %% UTILITY METHODS
    methods(Access=protected, Hidden)
        function error(obj, emsg, varargin)
            %ERROR Raise an error
            obj.isError = 1;
            if obj.batchMode
                error(emsg);
            else
                errordlg(emsg, varargin{:});
            end
        end

        function logOldP(obj, prm)
            %LOGOLDP Increment the old-style parameter counter
            %   collect stats on usage of old parameters
            if ~isKey(obj.oldPcount, prm)
                obj.oldPcount(prm) = 1;
            else
                obj.oldPcount(prm) = obj.oldPcount(prm) + 1;
            end
        end

        function loadParams(obj)
            obj.oldPcount = containers.Map();

            userParams = jrclust.utils.mToStruct(obj.configFile);
            if isfield(userParams, 'template_file') && ~isempty(userParams.template_file)
                try
                    t = jrclust.utils.mToStruct(jrclust.utils.absPath(userParams.template_file));
                    uParamNames = fieldnames(userParams);
                    for i = 1:numel(uParamNames) % merge s (specific) into t (general)
                        t.(uParamNames{i}) = userParams.(uParamNames{i});
                    end
                    userParams = t;
                catch ME
                    obj.warning(sprintf('Could not set template file %s: %s ', userParams.template_file, ME.message), 'Missing template file');
                end
            end

            uParamNames = fieldnames(userParams);
            unusedProps = cell(numel(uParamNames), 1);
            errorProps = cell(numel(uParamNames), 2);

            % loop through all user-defined parameters
            for i = 1:numel(uParamNames)
                % ignore configFile/template_file
                if ismember(uParamNames{i}, {'configFile', 'vcFile_prm', 'template_file'})
                    continue;
                elseif strcmp(uParamNames{i}, 'vcFile') % old single recording
                    if ~isempty(userParams.vcFile)
                        obj.setRawRecordings(userParams.vcFile);
                    end

                    continue;
                elseif strcmp(uParamNames{i}, 'csFile_merge')
                    if ~isempty(userParams.csFile_merge)
                        obj.setRawRecordings(userParams.csFile_merge);
                    end

                    continue;
                end

                % empty values in the param file take on their defaults
                if ~isempty(userParams.(uParamNames{i}))
                    if ~isprop(obj, uParamNames{i})
                        unusedProps{i} = uParamNames{i};
                        continue;
                    end

                    try
                        if ismember(uParamNames{i}, {'probe_file', 'probeFile'})
                            % try same directory as config file first
                            pf = jrclust.utils.absPath(userParams.(uParamNames{i}), fileparts(obj.configFile));
                            if isempty(pf) % fall back to default collection
                                pf = jrclust.utils.absPath(userParams.(uParamNames{i}), fullfile(jrclust.utils.basedir(), 'probes'));
                            end
                            obj.probeFile = pf;
                        end

                        obj.(uParamNames{i}) = userParams.(uParamNames{i});
                    catch ME % error or not a property we support
                        errorProps{i, 1} = uParamNames{i};
                        errorProps{i, 2} = ME.message;
                    end
                end
            end

            % warn user of unrecognized props
            uu = cellfun(@(u) ~isempty(u), unusedProps);
            if any(uu) && ~obj.batchMode
                wmsg = sprintf('The following properties were not recognized (possibly deprecated) and will be ignored:\n    %s', ...
                               strjoin(unusedProps(uu), '\n    '));
                obj.warning(wmsg, 'Unrecognized properties');
            end

            ee = cellfun(@(e) ~isempty(e), errorProps(:, 1));
            if any(ee) && ~obj.batchMode
                emsgs = arrayfun(@(i) strjoin(errorProps(i, :), ': '), find(ee), 'UniformOutput', 0);
                wmsg = sprintf('These properties were not set due to errors:\n* %s', strjoin(emsgs, '\n\n * '));
                obj.warning(wmsg, 'Unset properties');
            end

            if isempty(obj.rawRecordings)
                obj.error('Specify rawRecordings', 'No recording files');
                return;
            end

            % load probe from file
            if isempty(obj.probeFile)
                obj.error('Specify a probe file', 'Missing probe file');
                return;
            end
            try
                obj.loadProbe();
            catch ME
                obj.error(sprintf('%s', ME.message), 'Error loading probe file');
                return;
            end
        end

        % TODO: deprecate probe file and incorporate directly into config
        function loadProbe(obj)
            %LOADPROBE Load probe construct from .mat or .prb
            probe = doLoadProbe(obj.probeFile);
            fieldNames = fieldnames(probe);

            for i = 1:numel(fieldNames)
                fn = fieldNames{i};
                obj.(fn) = probe.(fn);
            end
        end

        function setRawRecordings(obj, rr)
            %SETRAWRECORDINGS Set rawRecordings field (depends on nonempty
            %   configFile)
            if ischar(rr)
                rr_ = jrclust.utils.absPath(rr, fileparts(obj.configFile));
                if ischar(rr_)
                    rr_ = {rr_};
                end
            elseif iscell(rr)
                rr_ = cellfun(@(rf) jrclust.utils.absPath(rf, fileparts(obj.configFile)), rr);
                rr_ = rr(~cellfun(@isempty, rr_));
            end

            if isempty(rr_)
                if iscell(rr)
                    rr = strjoin(rr, ', ');
                end

                obj.error(sprintf('Could not find recording(s): %s', rr));
                return;
            end

            % natsort filenames
            obj.rawRecordings = jrclust.utils.sortNat(rr_);
        end

        function validateParams(obj)
            %VALIDATEPARAMS Ensure parameters make sense
            % check probe is consistent
            if ~(jrclust.utils.ismatrixnum(obj.siteMap) && all(obj.siteMap > 0))
                obj.error('Malformed channel map', 'Bad probe configuration');
                return;
            end

            nc = numel(obj.siteMap);
            if ~(ismatrix(obj.siteLoc) && size(obj.siteLoc, 1) == nc)
                obj.error('Malformed probe geometry', 'Bad probe configuration');
                return;
            end

            if numel(obj.shankMap) ~= nc
                obj.error('Malformed shank indexing', 'Bad probe configuration');
                return;
            end

            obj.ignoreSites = obj.ignoreSites(ismember(obj.siteMap, obj.ignoreSites));

            % nSiteDir and/or nSitesExcl may not have been specified
            if isempty(obj.nSiteDir) || isempty(obj.nSitesExcl)
                siteDists = pdist2(obj.siteLoc, obj.siteLoc);

                % max over all sites of number of neighbors in merge radius
                nNeighMrg = max(sum(siteDists <= obj.evtMergeRad)); % 11/7/17 JJJ: med to max

                if isempty(obj.nSitesExcl)
                    % max over all sites of number of neighbors in detect radius
                    nNeighDetect = max(sum(siteDists <= obj.evtDetectRad)); % 11/7/17 JJJ: med to max
                    nsd = (nNeighDetect - 1)/2;
                    obj.nSitesExcl = nNeighDetect - nNeighMrg;
                else
                    nNeighDetect = nNeighMrg + obj.nSitesExcl;
                    nsd = (nNeighDetect - 1)/2;
                end

                if isempty(obj.nSiteDir)
                    obj.nSiteDir = nsd;
                end
            end

            if obj.nSitesEvt <= 0
                obj.error('nSitesExcl is too large or nSiteDir is too small', 'Bad configuration');
            end

            % try to infer a ground-truth file
            if isempty(obj.gtFile)
                try
                    obj.gtFile = jrclust.utils.subsExt(obj.configFile, '_gt.mat');
                catch % does not exist, leave empty
                end
            end

            obj.siteNeighbors = findSiteNeighbors(obj.siteLoc, 2*obj.nSiteDir + 1, obj.ignoreSites, obj.shankMap);

            % boost that gain
            obj.bitScaling = obj.bitScaling/obj.gainBoost;
        end

        function warning(obj, wmsg, varargin)
            %WARNING Raise a warning
            if obj.batchMode
                warning(wmsg);
            else
                warndlg(wmsg, varargin{:});
            end
        end
    end

    %% USER METHODS
    methods
        function edit(obj)
            %EDIT Edit the config file
            edit(obj.configFile);
        end

        function success = save(obj)
            %SAVE Write stored values to file
            success = 1;

            % first back up the old config file
            backupFile = jrclust.utils.subsExt(obj.configFile, '.prm.bak');
            try
                copyfile(obj.configFile, backupFile);
            catch ME % cowardly back out
                warning(ME.identifier, 'Could not back up old config file: %s', ME.message);
                success = 0;
                return;
            end

            try
                fid = fopen(obj.configFile, 'w');
            catch ME
                warning(ME.identifier, 'Could not open config file for writing: %s', ME.message);
                success = 0;
                return;
            end

            % TODO: read description from JSON file and organize by group
            fieldNames = sort(fieldnames(obj));
            for i = 1:numel(fieldNames)
                fn = fieldNames{i};
                % skip dependent properties as they can be computed on the fly
                p = findprop(obj, fn);
                if p.Dependent
                    continue;
                end
                fprintf(fid, '%s = %s;\n', fn, jrclust.utils.field2str(obj.(fn)));
            end

            fclose(fid);
        end

        function val = getOr(obj, fn, dv)
            %GETOR GET set value `obj.(fn)` OR default value `dv` if unset or empty
            if nargin < 3
                dv = [];
            end

            if ~isprop(obj, fn) || isempty(obj.(fn))
                val = dv;
            else
                val = obj.(fn);
            end
        end

        function success = makeCfg(obj, outputDir, binFiles, probeFile, templateFile, fAsk)
            %MAKECFG Make a config file
            success = 0;
            if obj.batchMode % no interactivity permitted
                return;
            end

            if nargin < 2
                outputDir = '';
            end
            if nargin < 3
                binFiles = {};
            end
            if nargin < 3
                probeFile = '';
            end
            if nargin < 4
                templateFile = '';
            end
            if nargin < 5
                fAsk = 1;
            end

            % set outputDir
            if isempty(outputDir)
                outputDir = uigetdir('', 'Select a directory to save your config file');
            end

            % set binFiles
            if ischar(binFiles)
                binFiles = {binFiles};
            end
            binFiles = cellfun(@jrclust.utils.absPath, binFiles, 'UniformOutput', 0);

            % absPath can return a cell if given a wildcard; handle mixed
            % output
            isWildcard = cellfun(@iscell, binFiles);
            if any(isWildcard)
                wildcards = binFiles(isWildcard);
                binFiles = binFiles(~isWildcard);
                for i = 1:numel(wildcards)
                    binFiles = [binFiles, wildcards{i}]; %#ok<AGROW>
                end
            end
            binFiles = unique(binFiles(~cellfun(@isempty, binFiles)));

            if isempty(binFiles)
                [fnames, fdir] = uigetfile({'*.bin'; '*.dat'}, 'Select a recording or recordings', 'MultiSelect', 'on');
                if ~ischar(fdir) % abort
                    return;
                end

                binFiles = fullfile(fdir, fnames);
                if ischar(binFiles)
                    binFiles = {binFiles};
                end
            end

            % TODO: deprecate probe file and incorporate directly into config
            probeFile = jrclust.utils.absPath(probeFile);
            if isempty(probeFile)
                [fname, fdir] = uigetfile({'*.prb'; '*.mat'}, ...
                                           'Select a probe file', ...
                                           fullfile(jrclust.utils.basedir(), 'probes'));
                if ~ischar(fdir) % abort
                    return;
                end

                probeFile = fullfile(fdir, fname);
            end

            % data file, template, and probe found; make the parameter file
            P = doCreateConfigFile(outputDir, binFiles, probeFile, templateFile, fAsk);
            if isempty(P)
                return;
            end

            fieldNames = intersect(fieldnames(obj),  fieldnames(P));
            for i = 1:numel(fieldNames)
                obj.(fieldNames{i}) = P.(fieldNames{i});
            end

            obj.flush(); % write to file
            success = 1;
        end

        function rd = recDurationSec(obj, recID)
            %RECDURATIONSECS Get duration of recording file(s) in seconds
            if nargin < 2 || isempty(recID)
                hRecs = cellfun(@(fn) jrclust.models.recording.Recording(fn, obj), obj.rawRecordings, 'UniformOutput', 0);
                rd = sum(cellfun(@(hR) hR.nSamples, hRecs))/obj.sampleRate;
            elseif recID < 1 || recID > numel(obj.rawRecordings)
                error('recording ID %d is invalid (there are %d recordings)', recID, numel(obj.rawRecordings));
            else
                hRec = jrclust.models.recording.Recording(obj.rawRecordings{recID}, obj);
                rd = hRec.nSamples/obj.sampleRate;
            end
        end

        function resetTemporaryParams(obj, prmKeys)
            %RESETTEMPORARYPARAMS Reset temporary parameters
            if nargin < 2 || isempty(prmKeys)
                prmKeys = keys(obj.tempParams);
            elseif nargin == 2
                if ischar(prmKeys)
                    prmKeys = {prmKeys};
                end
                % only try to reset parameters we actually have
                prmKeys = intersect(prmKeys, keys(obj.tempParams));
            end

            for i = 1:numel(prmKeys)
                fn = prmKeys{i};
                obj.(fn) = obj.tempParams(fn);
                remove(obj.tempParams, fn);
            end
        end

        function setTemporaryParams(obj, varargin)
            %SETTEMPORARYPARAMS Set temporary parameters to reset later
            prmKeys = varargin(1:2:end);
            prmVals = varargin(2:2:end);

            if numel(prmKeys) ~= numel(prmVals)
                warning('number of property names not equal to values; skipping');
                return;
            end

            for i = 1:numel(prmKeys)
                prmKey = prmKeys{i};

                % already set a temporary value for this parameter, reset
                % it or we'll lose the original
                if isKey(obj.tempParams, prmKey)
                    obj.resetTemporaryParams(prmKey);
                end
                try
                    obj.tempParams(prmKey) = obj.(prmKey); % save old value for later
                    obj.(prmKey) = prmVals{i};
                catch ME
                    remove(obj.tempParams, prmKey);
                    warning(ME.identifier, 'failed to set %s: %s', prmKey, ME.message);
                end
            end
        end
    end

    %% GETTERS/SETTERS
    methods
        % CARMode/vcCommonRef
        function set.CARMode(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            legalVals = {'mean', 'median', 'none'};
            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))
            obj.CARMode = val;
        end
        function val = get.vcCommonRef(obj)
            obj.logOldP('vcCommonRef');
            val = obj.CARMode;
        end
        function set.vcCommonRef(obj, val)
            obj.logOldP('vcCommonRef');
            obj.CARMode = val;
        end

        % autoMergeBy/autoMergeCriterion
        function set.autoMergeBy(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            legalVals = {'pearson', 'dist'};
            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))
            obj.autoMergeBy = val;
        end
        function val = get.autoMergeCriterion(obj)
            obj.logOldP('autoMergeCriterion');
            val = obj.autoMergeBy;
        end
        function set.autoMergeCriterion(obj, val)
            obj.logOldP('autoMergeCriterion');
            obj.autoMergeBy = val;
        end

        % auxChan/iChan_aux
        function set.auxChan(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.auxChan = val;
        end
        function val = get.iChan_aux(obj)
            obj.logOldP('iChan_aux');
            val = obj.auxChan;
        end
        function set.iChan_aux(obj, val)
            obj.logOldP('iChan_aux');
            obj.auxChan = val;
        end

        % auxFile/vcFile_aux
        function set.auxFile(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            hFun = @(x) ~isempty(jrclust.utils.absPath(x));
            assert(hFun(val));
            obj.auxFile = jrclust.utils.absPath(val);
        end
        function val = get.vcFile_aux(obj)
            obj.logOldP('vcFile_aux');
            val = obj.auxFile;
        end
        function set.vcFile_aux(obj, val)
            obj.logOldP('vcFile_aux');
            obj.auxFile = val;
        end

        % auxLabel/vcLabel_aux
        function set.auxLabel(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            obj.auxLabel = val;
        end
        function val = get.vcLabel_aux(obj)
            obj.logOldP('vcLabel_aux');
            val = obj.auxLabel;
        end
        function set.vcLabel_aux(obj, val)
            obj.logOldP('vcLabel_aux');
            obj.auxLabel = val;
        end

        % auxSampleRate/sRateHz_aux
        function set.auxSampleRate(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.auxSampleRate = val;
        end
        function val = get.sRateHz_aux(obj)
            obj.logOldP('sRateHz_aux');
            val = obj.auxSampleRate;
        end
        function set.sRateHz_aux(obj, val)
            obj.logOldP('sRateHz_aux');
            obj.auxSampleRate = val;
        end

        % auxScale/vrScale_aux
        function set.auxScale(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.auxScale = val;
        end
        function val = get.vrScale_aux(obj)
            obj.logOldP('vrScale_aux');
            val = obj.auxScale;
        end
        function set.vrScale_aux(obj, val)
            obj.logOldP('vrScale_aux');
            obj.auxScale = val;
        end

        % bitScaling/uV_per_bit
        function set.bitScaling(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.bitScaling = val;
        end
        function val = get.uV_per_bit(obj)
            obj.logOldP('uV_per_bit');
            val = obj.bitScaling;
        end
        function set.uV_per_bit(obj, val)
            obj.logOldP('uV_per_bit');
            obj.bitScaling = val;
        end

        % blankPeriod/blank_period_ms
        function set.blankPeriod(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.blankPeriod = val;
        end
        function val = get.blank_period_ms(obj)
            obj.logOldP('blank_period_ms');
            val = obj.blankPeriod;
        end
        function set.blank_period_ms(obj, val)
            obj.logOldP('blank_period_ms');
            obj.blankPeriod = val;
        end

        % blankThresh/blank_thresh
        function set.blankThresh(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.blankThresh = val;
        end
        function val = get.blank_thresh(obj)
            obj.logOldP('blank_thresh');
            val = obj.blankThresh;
        end
        function set.blank_thresh(obj, val)
            obj.logOldP('blank_thresh');
            obj.blankThresh = val;
        end

        % bytesPerSample
        function bp = get.bytesPerSample(obj)
            bp = jrclust.utils.typeBytes(obj.dataType);
        end

        % clusterFeature/vcFet
        function set.clusterFeature(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            if strcmp(val, 'gpca')
                val = 'pca';
                warning('gpca feature temporarily disabled; using pca');
            end
            legalVals = {'cov', 'energy', 'pca', 'vmin', 'vminmax', 'vpp'};
            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))
            obj.clusterFeature = val;
        end
        function val = get.vcFet(obj)
            obj.logOldP('vcFet');
            val = obj.clusterFeature;
        end
        function set.vcFet(obj, val)
            obj.logOldP('vcFet');
            obj.clusterFeature = val;
        end

        % colorMap/mrColor_proj
        function set.colorMap(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'nonnegative', 'numel', 9});
            hFun = @(x) reshape(x, 3, 3);
            val = hFun(val);
            obj.colorMap = val;
        end
        function val = get.mrColor_proj(obj)
            obj.logOldP('mrColor_proj');
            val = obj.colorMap;
        end
        function set.mrColor_proj(obj, val)
            obj.logOldP('mrColor_proj');
            obj.colorMap = val;
        end

        % configFile/vcFile_prm
        function set.configFile(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            hFun = @(x) ~isempty(jrclust.utils.absPath(x));
            assert(hFun(val));
            obj.configFile = jrclust.utils.absPath(val);
        end
        function cf = get.vcFile_prm(obj)
            obj.logOldP('vcFile_prm');
            cf = obj.configFile;
        end
        function set.vcFile_prm(obj, cf)
            obj.logOldP('vcFile_prm');
            obj.configFile = cf;
        end

        % corrRange/corrLim
        function set.corrRange(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'nonnegative', 'numel', 2, '<=', 1});
            obj.corrRange = val;
        end
        function val = get.corrLim(obj)
            obj.logOldP('corrLim');
            val = obj.corrRange;
        end
        function set.corrLim(obj, val)
            obj.logOldP('corrLim');
            obj.corrRange = val;
        end

        % dataType/vcDataType
        function set.dataType(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            legalVals = {'int16', 'uint16', 'int32', 'uint32', 'single', 'double'};
            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))
            obj.dataType = val;
        end
        function val = get.vcDataType(obj)
            obj.logOldP('vcDataType');
            val = obj.dataType;
        end
        function set.vcDataType(obj, val)
            obj.logOldP('vcDataType');
            obj.dataType = val;
        end

        % detectBipolar/fDetectBipolar
        function set.detectBipolar(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.detectBipolar = val;
        end
        function val = get.fDetectBipolar(obj)
            obj.logOldP('fDetectBipolar');
            val = obj.detectBipolar;
        end
        function set.fDetectBipolar(obj, val)
            obj.logOldP('fDetectBipolar');
            obj.detectBipolar = val;
        end

        % dispFeature/vcFet_show
        function set.dispFeature(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            legalVals = {'cov', 'pca', 'ppca', 'vpp'};
            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))
            obj.dispFeature = val;
        end
        function val = get.vcFet_show(obj)
            obj.logOldP('vcFet_show');
            val = obj.dispFeature;
        end
        function set.vcFet_show(obj, val)
            obj.logOldP('vcFet_show');
            obj.dispFeature = val;
        end

        % dispFilter/vcFilter_show
        function set.dispFilter(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            legalVals = {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'none'};
            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))
            obj.dispFilter = val;
        end
        function val = get.vcFilter_show(obj)
            obj.logOldP('vcFilter_show');
            val = obj.dispFilter;
        end
        function set.vcFilter_show(obj, val)
            obj.logOldP('vcFilter_show');
            obj.dispFilter = val;
        end

        % dispTimeLimits/tlim
        function set.dispTimeLimits(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'nonnegative', 'increasing', 'numel', 2});
            obj.dispTimeLimits = val;
        end
        function val = get.tlim(obj)
            obj.logOldP('tlim');
            val = obj.dispTimeLimits;
        end
        function set.tlim(obj, val)
            obj.logOldP('tlim');
            obj.dispTimeLimits = val;
        end

        % distCut/dc_percent
        function set.distCut(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'nonnegative', '<=', 100});
            obj.distCut = val;
        end
        function val = get.dc_percent(obj)
            obj.logOldP('dc_percent');
            val = obj.distCut;
        end
        function set.dc_percent(obj, val)
            obj.logOldP('dc_percent');
            obj.distCut = val;
        end

        % driftMerge/fDrift_merge
        function set.driftMerge(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.driftMerge = val;
        end
        function val = get.fDrift_merge(obj)
            obj.logOldP('fDrift_merge');
            val = obj.driftMerge;
        end
        function set.fDrift_merge(obj, val)
            obj.logOldP('fDrift_merge');
            obj.driftMerge = val;
        end

        % evtDetectRad/maxDist_site_spk_um
        function set.evtDetectRad(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.evtDetectRad = val;
        end
        function val = get.maxDist_site_spk_um(obj)
            obj.logOldP('maxDist_site_spk_um');
            val = obj.evtDetectRad;
        end
        function set.maxDist_site_spk_um(obj, val)
            obj.logOldP('maxDist_site_spk_um');
            obj.evtDetectRad = val;
        end

        % evtManualThreshSamp/spkThresh
        function val = get.evtManualThreshSamp(obj)
            val = obj.evtManualThresh / obj.bitScaling;
        end
        function val = get.spkThresh(obj)
            obj.logOldP('spkThresh');
            val = obj.evtManualThreshSamp;
        end

        % evtManualThresh/spkThresh_uV
        function set.evtManualThresh(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real'});
            obj.evtManualThresh = val;
        end
        function val = get.spkThresh_uV(obj)
            obj.logOldP('spkThresh_uV');
            val = obj.evtManualThresh;
        end
        function set.spkThresh_uV(obj, val)
            obj.logOldP('spkThresh_uV');
            obj.evtManualThresh = val;
        end

        % evtMergeRad/maxDist_site_um
        function set.evtMergeRad(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.evtMergeRad = val;
        end
        function val = get.maxDist_site_um(obj)
            obj.logOldP('maxDist_site_um');
            val = obj.evtMergeRad;
        end
        function set.maxDist_site_um(obj, val)
            obj.logOldP('maxDist_site_um');
            obj.evtMergeRad = val;
        end

        % evtWindow/spkLim_ms
        function set.evtWindow(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'increasing', 'numel', 2});
            hFun = @(x) diff(x) >= x(2);
            assert(hFun(val));
            obj.evtWindow = val;
        end
        function val = get.spkLim_ms(obj)
            obj.logOldP('spkLim_ms');
            val = obj.evtWindow;
        end
        function set.spkLim_ms(obj, val)
            obj.logOldP('spkLim_ms');
            obj.evtWindow = val;
        end

        % evtWindowMergeFactor/spkLim_factor_merge
        function set.evtWindowMergeFactor(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive', '<=', 1});
            obj.evtWindowMergeFactor = val;
        end
        function val = get.spkLim_factor_merge(obj)
            obj.logOldP('spkLim_factor_merge');
            val = obj.evtWindowMergeFactor;
        end
        function set.spkLim_factor_merge(obj, val)
            obj.logOldP('spkLim_factor_merge');
            obj.evtWindowMergeFactor = val;
        end

        % evtWindowRaw/spkLim_raw_ms
        function set.evtWindowRaw(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'increasing', 'numel', 2});
            hFun = @(x) diff(x) >= x(2);
            assert(hFun(val));
            obj.evtWindowRaw = val;
        end
        function val = get.spkLim_raw_ms(obj)
            obj.logOldP('spkLim_raw_ms');
            val = obj.evtWindowRaw;
        end
        function set.spkLim_raw_ms(obj, val)
            obj.logOldP('spkLim_raw_ms');
            obj.evtWindowRaw = val;
        end

        % evtWindowRawSamp/spkLim_raw
        function val = get.evtWindowRawSamp(obj)
            val = round(obj.evtWindowRaw * obj.sampleRate / 1000);
        end
        function set.evtWindowRawSamp(obj, val)
            obj.evtWindowRaw = val * 1000 / obj.sampleRate;
        end
        function val = get.spkLim_raw(obj)
            obj.logOldP('spkLim_raw');
            val = obj.evtWindowRawSamp;
        end
        function set.spkLim_raw(obj, val)
            obj.logOldP('spkLim_raw');
            obj.evtWindowRawSamp = val;
        end

        % evtWindowSamp/spkLim
        function val = get.evtWindowSamp(obj)
            val = round(obj.evtWindow * obj.sampleRate / 1000);
        end
        function set.evtWindowSamp(obj, val)
            obj.evtWindow = val * 1000 / obj.sampleRate;
        end
        function val = get.spkLim(obj)
            obj.logOldP('spkLim');
            val = obj.evtWindowSamp;
        end
        function set.spkLim(obj, val)
            obj.logOldP('spkLim');
            obj.evtWindowSamp = val;
        end

        % fftThresh/fft_thresh
        function set.fftThresh(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'scalar', 'nonnegative'});
            obj.fftThresh = val;
        end
        function val = get.fft_thresh(obj)
            obj.logOldP('fft_thresh');
            val = obj.fftThresh;
        end
        function set.fft_thresh(obj, val)
            obj.logOldP('fft_thresh');
            obj.fftThresh = val;
        end

        % filterType/vcFilter
        function set.filterType(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            legalVals = {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'none'};
            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))
            obj.filterType = val;
        end
        function val = get.vcFilter(obj)
            obj.logOldP('vcFilter');
            val = obj.filterType;
        end
        function set.vcFilter(obj, val)
            obj.logOldP('vcFilter');
            obj.filterType = val;
        end

        % filtOrder
        function set.filtOrder(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'scalar', 'positive'});
            obj.filtOrder = val;
        end

        % frFilterShape/filter_shape_rate
        function set.frFilterShape(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            legalVals = {'triangle', 'rectangle'};
            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))
            obj.frFilterShape = val;
        end
        function val = get.filter_shape_rate(obj)
            obj.logOldP('filter_shape_rate');
            val = obj.frFilterShape;
        end
        function set.filter_shape_rate(obj, val)
            obj.logOldP('filter_shape_rate');
            obj.frFilterShape = val;
        end

        % frPeriod/filter_sec_rate
        function set.frPeriod(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.frPeriod = val;
        end
        function val = get.filter_sec_rate(obj)
            obj.logOldP('filter_sec_rate');
            val = obj.frPeriod;
        end
        function set.filter_sec_rate(obj, val)
            obj.logOldP('filter_sec_rate');
            obj.frPeriod = val;
        end

        % frSampleRate/sRateHz_rate
        function set.frSampleRate(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.frSampleRate = val;
        end
        function val = get.sRateHz_rate(obj)
            obj.logOldP('sRateHz_rate');
            val = obj.frSampleRate;
        end
        function set.sRateHz_rate(obj, val)
            obj.logOldP('sRateHz_rate');
            obj.frSampleRate = val;
        end

        % freqLimBP/freqLim
        function set.freqLimBP(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'nonnegative', 'increasing', 'numel', 2});
            obj.freqLimBP = val;
        end
        function val = get.freqLim(obj)
            obj.logOldP('freqLim');
            val = obj.freqLimBP;
        end
        function set.freqLim(obj, val)
            obj.logOldP('freqLim');
            obj.freqLimBP = val;
        end

        % gainBoost/gain_boost
        function set.gainBoost(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.gainBoost = val;
        end
        function val = get.gain_boost(obj)
            obj.logOldP('gain_boost');
            val = obj.gainBoost;
        end
        function set.gain_boost(obj, val)
            obj.logOldP('gain_boost');
            obj.gainBoost = val;
        end

        % groupShank/fGroup_shank
        function set.groupShank(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.groupShank = val;
        end
        function val = get.fGroup_shank(obj)
            obj.logOldP('fGroup_shank');
            val = obj.groupShank;
        end
        function set.fGroup_shank(obj, val)
            obj.logOldP('fGroup_shank');
            obj.groupShank = val;
        end

        % gtFile/vcFile_gt
        function set.gtFile(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            hFun = @(x) ~isempty(jrclust.utils.absPath(x));
            assert(hFun(val));
            obj.gtFile = jrclust.utils.absPath(val);
        end
        function val = get.vcFile_gt(obj)
            obj.logOldP('vcFile_gt');
            val = obj.gtFile;
        end
        function set.vcFile_gt(obj, val)
            obj.logOldP('vcFile_gt');
            obj.gtFile = val;
        end

        % headerOffset/header_offset
        function set.headerOffset(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
            obj.headerOffset = val;
        end
        function val = get.header_offset(obj)
            obj.logOldP('header_offset');
            val = obj.headerOffset;
        end
        function set.header_offset(obj, val)
            obj.logOldP('header_offset');
            obj.headerOffset = val;
        end

        % ignoreSites/viSiteZero
        function set.ignoreSites(obj, val)
            validateattributes(val, {'numeric'}, {'integer', 'positive'});
            obj.ignoreSites = val;
        end
        function val = get.viSiteZero(obj)
            obj.logOldP('viSiteZero');
            val = obj.ignoreSites;
        end
        function set.viSiteZero(obj, val)
            obj.logOldP('viSiteZero');
            obj.ignoreSites = val;
        end

        % interpPC/fInterp_fet
        function set.interpPC(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.interpPC = val;
        end
        function val = get.fInterp_fet(obj)
            obj.logOldP('fInterp_fet');
            val = obj.interpPC;
        end
        function set.fInterp_fet(obj, val)
            obj.logOldP('fInterp_fet');
            obj.interpPC = val;
        end

        % lfpSampleRate/sRateHz_lfp
        function set.lfpSampleRate(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.lfpSampleRate = val;
        end
        function val = get.sRateHz_lfp(obj)
            obj.logOldP('sRateHz_lfp');
            val = obj.lfpSampleRate;
        end
        function set.sRateHz_lfp(obj, val)
            obj.logOldP('sRateHz_lfp');
            obj.lfpSampleRate = val;
        end

        % loadTimeLimits/tlim_load
        function set.loadTimeLimits(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'nonnegative', 'increasing', 'numel', 2});
            obj.loadTimeLimits = val;
        end
        function val = get.tlim_load(obj)
            obj.logOldP('tlim_load');
            val = obj.loadTimeLimits;
        end
        function set.tlim_load(obj, val)
            obj.logOldP('tlim_load');
            obj.loadTimeLimits = val;
        end

        % log10DeltaCut/delta1_cut
        function set.log10DeltaCut(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real'});
            obj.log10DeltaCut = val;
        end
        function val = get.delta1_cut(obj)
            obj.logOldP('delta1_cut');
            val = obj.log10DeltaCut;
        end
        function set.delta1_cut(obj, val)
            obj.logOldP('delta1_cut');
            obj.log10DeltaCut = val;
        end

        % log10RhoCut/rho_cut
        function set.log10RhoCut(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real'});
            obj.log10RhoCut = val;
        end
        function val = get.rho_cut(obj)
            obj.logOldP('rho_cut');
            val = obj.log10RhoCut;
        end
        function set.rho_cut(obj, val)
            obj.logOldP('rho_cut');
            obj.log10RhoCut = val;
        end

        % maxBytesLoad/MAX_BYTES_LOAD
        function set.maxBytesLoad(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.maxBytesLoad = val;
        end
        function val = get.MAX_BYTES_LOAD(obj)
            obj.logOldP('MAX_BYTES_LOAD');
            val = obj.maxBytesLoad;
        end
        function set.MAX_BYTES_LOAD(obj, val)
            obj.logOldP('MAX_BYTES_LOAD');
            obj.maxBytesLoad = val;
        end

        % maxClustersSite/maxCluPerSite
        function set.maxClustersSite(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.maxClustersSite = val;
        end
        function val = get.maxCluPerSite(obj)
            obj.logOldP('maxCluPerSite');
            val = obj.maxClustersSite;
        end
        function set.maxCluPerSite(obj, val)
            obj.logOldP('maxCluPerSite');
            obj.maxClustersSite = val;
        end

        % maxSecLoad/MAX_LOAD_SEC
        function set.maxSecLoad(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.maxSecLoad = val;
        end
        function val = get.MAX_LOAD_SEC(obj)
            obj.logOldP('MAX_LOAD_SEC');
            val = obj.maxSecLoad;
        end
        function set.MAX_LOAD_SEC(obj, val)
            obj.logOldP('MAX_LOAD_SEC');
            obj.maxSecLoad = val;
        end

        % maxUnitSim/maxWavCor
        function set.maxUnitSim(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'nonnegative', '<=', 1});
            obj.maxUnitSim = val;
        end
        function val = get.maxWavCor(obj)
            obj.logOldP('maxWavCor');
            val = obj.maxUnitSim;
        end
        function set.maxWavCor(obj, val)
            obj.logOldP('maxWavCor');
            obj.maxUnitSim = val;
        end

        % meanInterpFactor/nInterp_merge
        function set.meanInterpFactor(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.meanInterpFactor = val;
        end
        function val = get.nInterp_merge(obj)
            obj.logOldP('nInterp_merge');
            val = obj.meanInterpFactor;
        end
        function set.nInterp_merge(obj, val)
            obj.logOldP('nInterp_merge');
            obj.meanInterpFactor = val;
        end

        % minClusterSize/min_count
        function set.minClusterSize(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.minClusterSize = val;
        end
        function val = get.min_count(obj)
            obj.logOldP('min_count');
            val = obj.minClusterSize;
        end
        function set.min_count(obj, val)
            obj.logOldP('min_count');
            obj.minClusterSize = val;
        end

        % minNeighborsDetect/nneigh_min_detect
        function set.minNeighborsDetect(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'nonnegative', '<=', 2});
            obj.minNeighborsDetect = val;
        end
        function val = get.nneigh_min_detect(obj)
            obj.logOldP('nneigh_min_detect');
            val = obj.minNeighborsDetect;
        end
        function set.nneigh_min_detect(obj, val)
            obj.logOldP('nneigh_min_detect');
            obj.minNeighborsDetect = val;
        end

        % minSitesWeightFeatures/min_sites_mask
        function set.minSitesWeightFeatures(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.minSitesWeightFeatures = val;
        end
        function val = get.min_sites_mask(obj)
            obj.logOldP('min_sites_mask');
            val = obj.minSitesWeightFeatures;
        end
        function set.min_sites_mask(obj, val)
            obj.logOldP('min_sites_mask');
            obj.minSitesWeightFeatures = val;
        end

        % nClusterIntervals/nTime_clu
        function set.nClusterIntervals(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nClusterIntervals = val;
        end
        function val = get.nTime_clu(obj)
            obj.logOldP('nTime_clu');
            val = obj.nClusterIntervals;
        end
        function set.nTime_clu(obj, val)
            obj.logOldP('nTime_clu');
            obj.nClusterIntervals = val;
        end

        % nClustersShowAux/nClu_show_aux
        function set.nClustersShowAux(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nClustersShowAux = val;
        end
        function val = get.nClu_show_aux(obj)
            obj.logOldP('nClu_show_aux');
            val = obj.nClustersShowAux;
        end
        function set.nClu_show_aux(obj, val)
            obj.logOldP('nClu_show_aux');
            obj.nClustersShowAux = val;
        end

        % nDiffOrder/nDiff_filt
        function set.nDiffOrder(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'nonnegative', '<=', 3});
            obj.nDiffOrder = val;
        end
        function val = get.nDiff_filt(obj)
            obj.logOldP('nDiff_filt');
            val = obj.nDiffOrder;
        end
        function set.nDiff_filt(obj, val)
            obj.logOldP('nDiff_filt');
            obj.nDiffOrder = val;
        end

        % nLoadsMaxPreview/nLoads_max_preview
        function set.nLoadsMaxPreview(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nLoadsMaxPreview = val;
        end
        function val = get.nLoads_max_preview(obj)
            obj.logOldP('nLoads_max_preview');
            val = obj.nLoadsMaxPreview;
        end
        function set.nLoads_max_preview(obj, val)
            obj.logOldP('nLoads_max_preview');
            obj.nLoadsMaxPreview = val;
        end

        % nPCsPerSite/nPcPerChan
        function set.nPCsPerSite(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive', '<=', 3});
            obj.nPCsPerSite = val;
        end
        function val = get.nPcPerChan(obj)
            obj.logOldP('nPcPerChan');
            val = obj.nPCsPerSite;
        end
        function set.nPcPerChan(obj, val)
            obj.logOldP('nPcPerChan');
            obj.nPCsPerSite = val;
        end

        % nPassesMerge/nRepeat_merge
        function set.nPassesMerge(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nPassesMerge = val;
        end
        function val = get.nRepeat_merge(obj)
            obj.logOldP('nRepeat_merge');
            val = obj.nPassesMerge;
        end
        function set.nRepeat_merge(obj, val)
            obj.logOldP('nRepeat_merge');
            obj.nPassesMerge = val;
        end

        % nPeaksFeatures/nFet_use
        function set.nPeaksFeatures(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', '<=', 3});
            obj.nPeaksFeatures = val;
        end
        function val = get.nFet_use(obj)
            obj.logOldP('nFet_use');
            val = obj.nPeaksFeatures;
        end
        function set.nFet_use(obj, val)
            obj.logOldP('nFet_use');
            obj.nPeaksFeatures = val;
        end

        % nSamplesPad/nPad_filt
        function set.nSamplesPad(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nSamplesPad = val;
        end
        function val = get.nPad_filt(obj)
            obj.logOldP('nPad_filt');
            val = obj.nSamplesPad;
        end
        function set.nPad_filt(obj, val)
            obj.logOldP('nPad_filt');
            obj.nSamplesPad = val;
        end

        % nSecsLoadPreview/sec_per_load_preview
        function set.nSecsLoadPreview(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.nSecsLoadPreview = val;
        end
        function val = get.sec_per_load_preview(obj)
            obj.logOldP('sec_per_load_preview');
            val = obj.nSecsLoadPreview;
        end
        function set.sec_per_load_preview(obj, val)
            obj.logOldP('sec_per_load_preview');
            obj.nSecsLoadPreview = val;
        end

        % nSegmentsTraces/nTime_traces
        function set.nSegmentsTraces(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nSegmentsTraces = val;
        end
        function val = get.nTime_traces(obj)
            obj.logOldP('nTime_traces');
            val = obj.nSegmentsTraces;
        end
        function set.nTime_traces(obj, val)
            obj.logOldP('nTime_traces');
            obj.nSegmentsTraces = val;
        end

        % nSiteDir/maxSite
        function set.nSiteDir(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nSiteDir = val;
        end
        function val = get.maxSite(obj)
            obj.logOldP('maxSite');
            val = obj.nSiteDir;
        end
        function set.maxSite(obj, val)
            obj.logOldP('maxSite');
            obj.nSiteDir = val;
        end

        % nSites
        function ns = get.nSites(obj)
            ns = numel(obj.siteMap);
        end

        % nSitesEvt
        function ns = get.nSitesEvt(obj)
            ns = 2*obj.nSiteDir - obj.nSitesExcl + 1;
        end

        % nSitesExcl/nSites_ref
        function set.nSitesExcl(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nSitesExcl = val;
        end
        function val = get.nSites_ref(obj)
            obj.logOldP('nSites_ref');
            val = obj.nSitesExcl;
        end
        function set.nSites_ref(obj, val)
            obj.logOldP('nSites_ref');
            obj.nSitesExcl = val;
        end

        % nSkip/nSkip_show
        function set.nSkip(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nSkip = val;
        end
        function val = get.nSkip_show(obj)
            obj.logOldP('nSkip_show');
            val = obj.nSkip;
        end
        function set.nSkip_show(obj, val)
            obj.logOldP('nSkip_show');
            obj.nSkip = val;
        end

        % nSpikesFigProj/nShow_proj
        function set.nSpikesFigProj(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nSpikesFigProj = val;
        end
        function val = get.nShow_proj(obj)
            obj.logOldP('nShow_proj');
            val = obj.nSpikesFigProj;
        end
        function set.nShow_proj(obj, val)
            obj.logOldP('nShow_proj');
            obj.nSpikesFigProj = val;
        end

        % nSpikesFigWav/nSpk_show
        function set.nSpikesFigWav(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nSpikesFigWav = val;
        end
        function val = get.nSpk_show(obj)
            obj.logOldP('nSpk_show');
            val = obj.nSpikesFigWav;
        end
        function set.nSpk_show(obj, val)
            obj.logOldP('nSpk_show');
            obj.nSpikesFigWav = val;
        end

        % nThreadsGPU/nThreads
        function set.nThreadsGPU(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.nThreadsGPU = val;
        end
        function val = get.nThreads(obj)
            obj.logOldP('nThreads');
            val = obj.nThreadsGPU;
        end
        function set.nThreads(obj, val)
            obj.logOldP('nThreads');
            obj.nThreadsGPU = val;
        end

        % outlierThresh/thresh_mad_clu
        function set.outlierThresh(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.outlierThresh = val;
        end
        function val = get.thresh_mad_clu(obj)
            obj.logOldP('thresh_mad_clu');
            val = obj.outlierThresh;
        end
        function set.thresh_mad_clu(obj, val)
            obj.logOldP('thresh_mad_clu');
            obj.outlierThresh = val;
        end

        % outputDir
        function val = get.outputDir(obj)
            if isempty(obj.outputDir)
                val = fileparts(obj.configFile);
            else
                val = obj.outputDir;
            end
        end
        function set.outputDir(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            hFun = @(x) ~isempty(jrclust.utils.absPath(x));
            assert(hFun(val));
            obj.outputDir = jrclust.utils.absPath(val);
        end

        % probeFile/probe_file
        function set.probeFile(obj, val)
            % first check where the config file is located
            val_ = jrclust.utils.absPath(val);
            % if we can't find it there, try the standard location
            if isempty(val_)
                val_ = jrclust.utils.absPath(val, fullfile(jrclust.utils.basedir(), 'probes'));
            end
            assert(isfile(val_), 'could not find probe file ''%s''', val);
            obj.probeFile = val_;
        end
        function val = get.probe_file(obj)
            obj.logOldP('probe_file');
            val = obj.probeFile;
        end
        function set.probe_file(obj, val)
            obj.logOldP('probe_file');
            obj.probeFile = val;
        end

        % probePad/vrSiteHW
        function set.probePad(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'positive', 'numel', 2});
            obj.probePad = val;
        end
        function val = get.vrSiteHW(obj)
            obj.logOldP('vrSiteHW');
            val = obj.probePad;
        end
        function set.vrSiteHW(obj, val)
            obj.logOldP('vrSiteHW');
            obj.probePad = val;
        end

        % projTimeLimits/tLimFigProj
        function set.projTimeLimits(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'nonnegative', 'increasing', 'numel', 2});
            obj.projTimeLimits = val;
        end
        function val = get.tLimFigProj(obj)
            obj.logOldP('tLimFigProj');
            val = obj.projTimeLimits;
        end
        function set.tLimFigProj(obj, val)
            obj.logOldP('tLimFigProj');
            obj.projTimeLimits = val;
        end

        % psthTimeBin/tbin_psth
        function set.psthTimeBin(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.psthTimeBin = val;
        end
        function val = get.tbin_psth(obj)
            obj.logOldP('tbin_psth');
            val = obj.psthTimeBin;
        end
        function set.tbin_psth(obj, val)
            obj.logOldP('tbin_psth');
            obj.psthTimeBin = val;
        end

        % psthTimeLimits/tlim_psth
        function set.psthTimeLimits(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'increasing', 'numel', 2});
            hFun = @(x) diff(x) >= x(2);
            assert(hFun(val));
            obj.psthTimeLimits = val;
        end
        function val = get.tlim_psth(obj)
            obj.logOldP('tlim_psth');
            val = obj.psthTimeLimits;
        end
        function set.tlim_psth(obj, val)
            obj.logOldP('tlim_psth');
            obj.psthTimeLimits = val;
        end

        % psthXTick/xtick_psth
        function set.psthXTick(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.psthXTick = val;
        end
        function val = get.xtick_psth(obj)
            obj.logOldP('xtick_psth');
            val = obj.psthXTick;
        end
        function set.xtick_psth(obj, val)
            obj.logOldP('xtick_psth');
            obj.psthXTick = val;
        end

        % ramToGPUFactor/nLoads_gpu
        function set.ramToGPUFactor(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.ramToGPUFactor = val;
        end
        function val = get.nLoads_gpu(obj)
            obj.logOldP('nLoads_gpu');
            val = obj.ramToGPUFactor;
        end
        function set.nLoads_gpu(obj, val)
            obj.logOldP('nLoads_gpu');
            obj.ramToGPUFactor = val;
        end

        % randomSeed
        function set.randomSeed(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
            obj.randomSeed = val;
        end

        % refracInt/spkRefrac_ms
        function set.refracInt(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.refracInt = val;
        end
        function val = get.spkRefrac_ms(obj)
            obj.logOldP('spkRefrac_ms');
            val = obj.refracInt;
        end
        function set.spkRefrac_ms(obj, val)
            obj.logOldP('spkRefrac_ms');
            obj.refracInt = val;
        end

        % refracIntSamp/spkRefrac
        function val = get.refracIntSamp(obj)
            val = round(obj.refracInt * obj.sampleRate / 1000);
        end
        function set.refracIntSamp(obj, val)
            obj.refracInt = val * 1000 / obj.sampleRate;
        end
        function val = get.spkRefrac(obj)
            obj.logOldP('spkRefrac');
            val = obj.refracIntSamp;
        end
        function set.spkRefrac(obj, val)
            obj.logOldP('spkRefrac');
            obj.refracIntSamp = val;
        end

        % sampleRate/sRateHz
        function set.sampleRate(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'});
            obj.sampleRate = val;
        end
        function val = get.sRateHz(obj)
            obj.logOldP('sRateHz');
            val = obj.sampleRate;
        end
        function set.sRateHz(obj, val)
            obj.logOldP('sRateHz');
            obj.sampleRate = val;
        end

        % sessionName
        function val = get.sessionName(obj)
            [~, val, ~] = fileparts(obj.configFile);
        end

        % shankMap/viShank_site
        function set.shankMap(obj, val)
            validateattributes(val, {'numeric'}, {'integer', 'positive'});
            obj.shankMap = val;
        end
        function val = get.viShank_site(obj)
            obj.logOldP('viShank_site');
            val = obj.shankMap;
        end
        function set.viShank_site(obj, val)
            obj.logOldP('viShank_site');
            obj.shankMap = val;
        end

        % showRaw/fWav_raw_show
        function set.showRaw(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.showRaw = val;
        end
        function val = get.fWav_raw_show(obj)
            obj.logOldP('fWav_raw_show');
            val = obj.showRaw;
        end
        function set.fWav_raw_show(obj, val)
            obj.logOldP('fWav_raw_show');
            obj.showRaw = val;
        end

        % showSpikeCount/fText
        function set.showSpikeCount(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.showSpikeCount = val;
        end
        function val = get.fText(obj)
            obj.logOldP('fText');
            val = obj.showSpikeCount;
        end
        function set.fText(obj, val)
            obj.logOldP('fText');
            obj.showSpikeCount = val;
        end

        % siteCorrThresh/thresh_corr_bad_site
        function set.siteCorrThresh(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'nonnegative', '<', 1});
            obj.siteCorrThresh = val;
        end
        function val = get.thresh_corr_bad_site(obj)
            obj.logOldP('thresh_corr_bad_site');
            val = obj.siteCorrThresh;
        end
        function set.thresh_corr_bad_site(obj, val)
            obj.logOldP('thresh_corr_bad_site');
            obj.siteCorrThresh = val;
        end

        % siteLoc/mrSiteXY
        function set.siteLoc(obj, val)
            validateattributes(val, {'numeric'}, {'real', 'nonnegative'});
            obj.siteLoc = val;
        end
        function val = get.mrSiteXY(obj)
            obj.logOldP('mrSiteXY');
            val = obj.siteLoc;
        end
        function set.mrSiteXY(obj, val)
            obj.logOldP('mrSiteXY');
            obj.siteLoc = val;
        end

        % siteMap/viSite2Chan
        function set.siteMap(obj, val)
            validateattributes(val, {'numeric'}, {'integer', 'positive'});
            obj.siteMap = val;
        end
        function val = get.viSite2Chan(obj)
            obj.logOldP('viSite2Chan');
            val = obj.siteMap;
        end
        function set.viSite2Chan(obj, val)
            obj.logOldP('viSite2Chan');
            obj.siteMap = val;
        end

        % siteNeighbors/miSites
        function set.siteNeighbors(obj, val) % danger zone: don't set this manually
            validateattributes(val, {'numeric'}, {'integer', 'positive'});
            obj.siteNeighbors = val;
        end
        function val = get.miSites(obj)
            obj.logOldP('miSites');
            val = obj.siteNeighbors;
        end
        function set.miSites(obj, val)
            obj.logOldP('miSites');
            obj.siteNeighbors = val;
        end

        % spikeThreshMax/spkThresh_max_uV
        function set.spikeThreshMax(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'});
            obj.spikeThreshMax = val;
        end
        function val = get.spkThresh_max_uV(obj)
            obj.logOldP('spkThresh_max_uV');
            val = obj.spikeThreshMax;
        end
        function set.spkThresh_max_uV(obj, val)
            obj.logOldP('spkThresh_max_uV');
            obj.spikeThreshMax = val;
        end

        % tallSkinny/fTranspose_bin
        function set.tallSkinny(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.tallSkinny = val;
        end
        function val = get.fTranspose_bin(obj)
            obj.logOldP('fTranspose_bin');
            val = obj.tallSkinny;
        end
        function set.fTranspose_bin(obj, val)
            obj.logOldP('fTranspose_bin');
            obj.tallSkinny = val;
        end

        % threshFile/vcFile_thresh
        function set.threshFile(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            hFun = @(x) ~isempty(jrclust.utils.absPath(x));
            assert(hFun(val));
            obj.threshFile = jrclust.utils.absPath(val);
        end
        function val = get.vcFile_thresh(obj)
            obj.logOldP('vcFile_thresh');
            val = obj.threshFile;
        end
        function set.vcFile_thresh(obj, val)
            obj.logOldP('vcFile_thresh');
            obj.threshFile = val;
        end

        % trialFile/vcFile_trial
        function set.trialFile(obj, val)
            validateattributes(val, {'char'}, {'scalartext'});
            hFun = @(x) ~isempty(jrclust.utils.absPath(x));
            assert(hFun(val));
            obj.trialFile =jrclust.utils.absPath(val);
        end
        function val = get.vcFile_trial(obj)
            obj.logOldP('vcFile_trial');
            val = obj.trialFile;
        end
        function set.vcFile_trial(obj, val)
            obj.logOldP('vcFile_trial');
            obj.trialFile = val;
        end

        % umPerPix/um_per_pix
        function set.umPerPix(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'positive', 'integer'});
            obj.umPerPix = val;
        end
        function val = get.um_per_pix(obj)
            obj.logOldP('um_per_pix');
            val = obj.umPerPix;
        end
        function set.um_per_pix(obj, val)
            obj.logOldP('um_per_pix');
            obj.umPerPix = val;
        end

        % useElliptic/fEllip
        function set.useElliptic(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.useElliptic = val;
        end
        function val = get.fEllip(obj)
            obj.logOldP('fEllip');
            val = obj.useElliptic;
        end
        function set.fEllip(obj, val)
            obj.logOldP('fEllip');
            obj.useElliptic = val;
        end

        % useGPU/fGpu
        function set.useGPU(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.useGPU = val;
        end
        function val = get.fGpu(obj)
            obj.logOldP('fGpu');
            val = obj.useGPU;
        end
        function set.fGpu(obj, val)
            obj.logOldP('fGpu');
            obj.useGPU = val;
        end

        % useGlobalDistCut/fDc_global
        function set.useGlobalDistCut(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.useGlobalDistCut = val;
        end
        function val = get.fDc_global(obj)
            obj.logOldP('fDc_global');
            val = obj.useGlobalDistCut;
        end
        function set.fDc_global(obj, val)
            obj.logOldP('fDc_global');
            obj.useGlobalDistCut = val;
        end

        % useParfor/fParfor
        function set.useParfor(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.useParfor = val;
        end
        function val = get.fParfor(obj)
            obj.logOldP('fParfor');
            val = obj.useParfor;
        end
        function set.fParfor(obj, val)
            obj.logOldP('fParfor');
            obj.useParfor = val;
        end

        % userFiltKernel/vnFilter_user
        function set.userFiltKernel(obj, val)
            validateattributes(val, {'numeric'}, {'integer'});
            obj.userFiltKernel = val;
        end
        function val = get.vnFilter_user(obj)
            obj.logOldP('vnFilter_user');
            val = obj.userFiltKernel;
        end
        function set.vnFilter_user(obj, val)
            obj.logOldP('vnFilter_user');
            obj.userFiltKernel = val;
        end

        % verbose/fVerbose
        function set.verbose(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.verbose = val;
        end
        function val = get.fVerbose(obj)
            obj.logOldP('fVerbose');
            val = obj.verbose;
        end
        function set.fVerbose(obj, val)
            obj.logOldP('fVerbose');
            obj.verbose = val;
        end

        % weightFeatures/fSpatialMask_clu
        function set.weightFeatures(obj, val)
            validateattributes(val, {'logical', 'double'}, {'scalar'});
            hFun = @(x) logical(x);
            val = hFun(val);
            obj.weightFeatures = val;
        end
        function val = get.fSpatialMask_clu(obj)
            obj.logOldP('fSpatialMask_clu');
            val = obj.weightFeatures;
        end
        function set.fSpatialMask_clu(obj, val)
            obj.logOldP('fSpatialMask_clu');
            obj.weightFeatures = val;
        end
    end
end
