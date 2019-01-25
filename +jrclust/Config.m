classdef Config < dynamicprops
    %CONFIG JRCLUST session configuration
    % replacement for P struct

    %% OBJECT-LEVEL PROPERTIES
    properties (Hidden, SetAccess=private, SetObservable)
        errMsg;
        isLoaded;
        isError;
        oldPcount;
        paramSet;
        oldParamSet;
        tempParams;
    end

    %% OLD-STYLE PARAMS, publicly settable (will be deprecated after a grace period)
    properties (SetObservable, Dependent, Hidden, Transient)
        autoMergeCriterion;         % => autoMergeBy
        blank_thresh;               % => blankThresh
        csFile_merge;               % => rawRecordings
        delta1_cut;                 % => log10DeltaCut
        fft_thresh;                 % => fftThreshMad
        fGpu;                       % => useGPU
        filter_sec_rate;            % => frPeriod
        filter_shape_rate;          % => frFilterShape
        fVerbose;                   % => verbose
        fWav_raw_show;              % => showRaw
        gain_boost;                 % => gainBoost
        header_offset;              % => headerOffset
        iChan_aux;                  % => auxChan
        MAX_BYTES_LOAD;             % => maxBytesLoad
        MAX_LOAD_SEC;               % => maxSecLoad
        maxCluPerSite;              % => maxClustersSite
        maxDist_site_um;            % => evtMergeRad
        maxDist_site_spk_um;        % => evtDetectRad
        maxSite;                    % => nSiteDir
        min_count;                  % => minClusterSize
        mrSiteXY;                   % => siteLoc
        nClu_show_aux;              % => nClustersShowAux
        nLoads_gpu;                 % => ramToGPUFactor
        nPad_filt;                  % => nSamplesPad
        nRepeat_merge;              % => nPassesMerge
        nSites_ref;                 % => nSitesExcl
        nSkip_lfp;                  % => lfpDsFactor
        probe_file;                 % => probeFile
        rho_cut;                    % => log10RhoCut
        spkLim;                     % => evtWindowSamp
        spkLim_ms;                  % => evtWindow
        spkLim_raw;                 % => evtWindowRawSamp
        spkLim_raw_ms;              % => evtWindowRaw
        spkRefrac;                  % => refracIntSamp
        spkRefrac_ms;               % => refracInt
        spkThresh;                  % => evtManualThreshSamp
        spkThresh_uV;               % => evtManualThresh
        sRateHz;                    % => sampleRate
        sRateHz_aux;                % => auxSampleRate
        sRateHz_lfp;                % => lfpSampleRate
        sRateHz_rate;               % => frSampleRate
        thresh_corr_bad_site;       % => siteCorrThresh
        thresh_mad_clu;             % => outlierThresh
        tlim;                       % => dispTimeLimits
        tlim_load;                  % => loadTimeLimits
        uV_per_bit;                 % => bitScaling
        vcCommonRef;                % => CARMode
        vcDataType;                 % => dataType
        vcDetrend_postclu;          % => RDDetrendMode
        vcFet;                      % => clusterFeature
        vcFet_show;                 % => dispFeature
        vcFile_aux;                 % => auxFile
        vcFile_gt;                  % => gtFile
        vcFile_prm;                 % => configFile
        vcFile_thresh;              % => threshFile
        vcFile_trial;               % => trialFile
        vcFilter;                   % => filterType
        vcFilter_show;              % => dispFilter
        vcLabel_aux;                % => auxLabel
        viShank_site;               % => shankMap
        vnFilter_user;              % => userFiltKernel
        vrSiteHW;                   % => probePad
        viSite2Chan;                % => siteMap
        viSiteZero;                 % => ignoreSites
        vrScale_aux;                % => auxScale
    end

    %% OLD-STLYE PARAMS, not publicly settable
    properties (Dependent, SetAccess=private, Hidden)
        miSites;                    % => siteNeighbors
    end

    %% NEW-STYLE PARAMS, publicly settable
    properties (SetObservable)
        batchMode = 0;              % suppress *all* messages if true
        verbose = 1;                % be chatty while processing

        % computation params
        useParfor = 1;                % use parfor where appropriate
        gpuLoadFactor = 5;          % GPU memory usage factor (4x means 1/4 of GPU memory can be loaded)
        nThreadsGPU = 128;             % number of gpu threads
        randomSeed = 0;             % random seed
        ramToGPUFactor = 8;         % ratio: RAM / (GPU memory) (increase this number if GPU memory error)
        useGPU = 1;                 % use GPU in computation if true

        % file location params
        outputDir = '';             % directory in which to place output files

        % recording params
        bitScaling = 0.30518;       % bit scaling factor (uV/bit)
        configFile;                 % parameter file
        dataType = 'int16';            % raw data binary format
        gainBoost = 1;              % multiply the raw recording by this gain to boost uV/bit
        gtFile = '';                % ground truth file (default: SESSION_NAME_gt.mat) (TODO: specify format)
        headerOffset = 0;           % file header offset, in bytes
        lfpSampleRate = 2500;       % sample rate of the LFP recording, in Hz
        nChans = 385;               % number of channels stored in recording
        probeFile;                  % probe file to use (.prb, .mat)
        probePad;                   %
        rawRecordings;              % collection of recordings
        sampleRate = 30000;         % sample rate of the recording, in Hz
        shankMap;                   % index of shank to which a site belongs
        siteLoc;                    % x-y locations of channels on the probe, in microns
        siteMap;                    % channel mapping; row i in the data corresponds to channel `siteMap(i)`

        % preprocessing params
        useElliptic = 1;            % use elliptic filter if 1 (and only if filterType='bandpass')
        fftThresh = 0;              % automatically remove frequency outliers (unit:MAD, 10 recommended, 0 to disable). Verify by running "jrc traces" and press "p" to view the power spectrum.
        filtOrder = 3;              % bandpass filter order
        filterType = 'ndiff';       % filter to use {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'fftdiff', 'none'}
        freqLimBP = [300 3000];       % frequency cut-off limit for filterType='bandpass' (ignored otherwise)
        freqLimNotch = [];
        freqLimStop = [];
        fTranspose_bin = 1;
        loadTimeLimits = [];        % time range of recording to load, in s (use whole range if empty)
        maxBytesLoad = [];          % default memory loading block size (bytes)
        maxSecLoad = [];            % maximum loading duration (seconds) (overrides 'maxBytesLoad')
        nSamplesPad = 100;          % number of samples to overlap between multiple loading (filter edge safe)
        userFiltKernel = [];        % custom filter kernel (optional unless filterType='user')
        CARMode = 'mean';           % common average referencing mode (one of 'none', 'mean', 'median')

        % spike detection params
        blankPeriod = 5;        % (miliseconds) Duration of blanking when the common mean exceeds a threhold (blank_thresh)
        blankThresh = [];           % reject spikes exceeding the channel mean after filtering (MAD unit), ignored if [] or 0
        evtDetectRad = 75;          % radius for extracting waveforms, in microns (used if nSiteDir and nSitesExcl are empty)
        evtManualThresh = [];     % manual spike detection threshold, in microvolts
        evtMergeRad = 50;           % radius of spike event merging, in microns
        evtWindow = [-0.25 0.75];   % interval around event to extract filtered spike waveforms, in ms
        evtWindowRaw = [-0.5 1.5];  % interval around event to extract raw spike waveforms, in ms
        groupShank = 0;           % group all sites in the same shank if true
        ignoreSites = [];           % sites to manually ignore in the sorting
        nDiffOrder = 2;             % Differentiation filter for filterType='sgdiff', ignored otherwise. Set to [] to disable. 2n+1 samples used for centered differentiation
        minNeighborsDetect = 0;      % Min. number of neighbors near the spike below threshold. choose between [0,1,2]
        nSiteDir;                   % number of neighboring sites to group in each direction (TODO: deprecate this)
        nSitesExcl;                 % number of sites to exclude from the spike waveform group
        qqFactor = 5;
        refracInt = 0.25;         % spike refractory interval, in ms
        siteCorrThresh = 0;         % reject bad sites based on max correlation with neighboring sites, using raw waveforms; ignored if 0
        spkThresh_max_uV = [];      % maximum absolute amp. allowed
        threshFile = '';            % path to .mat file storing spike detection thresholds (created by 'preview' GUI)

        % feature extraction params
        clusterFeature = 'pca';             % feature to use in clustering
        interpPC = 1;                    % interpolate waveforms for feature projection to find optimal delay (2x interp) if true
        fSpatialMask_clu = 0;               % apply spatial mask calculated from the distances between sites to the peak site (half-scale: evtDetectRad)
        minSitesWeightFeatures = 5;                 % minimum number of sites to have to apply spatial mask
        nPeaksFeatures = 2;                       % undocumented
        nPCsPerSite = 1;
        timeFeatureFactor;                % undocumented

        % clustering params
        autoMergeBy = 'pearson';            % metric to use when automerging clusters
        distCut = 2;                        % percentile at which to cut off distance in rho computation
        driftMerge = 1;                   % compute multiple waveforms at three drift locations based on the spike position if true
        log10DeltaCut = 0.6;                % the base-10 log of the delta cutoff value
        log10RhoCut = -2.5;                 % the base-10 log of the rho cutoff value
        maxClustersSite = 20;               % maximum number of clusters per site if local detrending is used
        minClusterSize = 30;                % minimum cluster size (set to 2*#features if lower)
        maxUnitSim = 0.98;                   %
        nInterp_merge = 1;                  % Interpolation factor for the mean unit waveforms, set to 1 to disable
        nPassesMerge = 10;                  % number of passes for unit mean raw waveform-based merging
        outlierThresh = 7.5;                % threshold to remove outlier spikes for each cluster, in MAD
        nClusterIntervals = 4;                      % number of time periods over which to cluster separately (later to be merged after clustering)
        RDDetrendMode = 'global';           %
        evtWindowMergeFactor = 1;            % Waveform range for computing the correlation. evtWindowMergeFactor <= spkLim_raw_factor_merge. circa v3.1.8
        useGlobalDistCut = 0;               % use a global distance cutoff for all sites if true; otherwise use different cutoff values for each site

        % display params
        corrRange = [0.9 1];
        dispFeature = 'vpp';                % feature to display in time/projection views
        dispFilter = '';                    %
        dispTimeLimits = [0 0.2];           % time range to display (in seconds)
        figList;                            % figure handles to display in manual view
        showSpikeCount = 1;                          %
        maxAmp = 250;
        colorMap = [213 219 235; ...
                    0   130 196; ...
                    240 119 22]/256;
        nShow = 200;                        % maximum number of traces to show [D?# spikes to show]
        nSpikesFigProj = 500;                   % maximum number of features to show in projection
        nSitesFigProj = 5;                  % number of sites to display in the feature projection view
        nSkip = 1;
        nSegmentsTraces = 1;                   % number of time segments to display. Set to 1 to show one continuous time segment
        nSpikesFigWav = 30;                     % show spike waveforms for manual clustering
        pcPair = [1 2];                     % PC projection to show (1 vs 2; 1 vs 3; 2 vs 3), can be toggled
        showRaw = 0;                        % show raw waveforms in main view if true
        projTimeLimits = [];                   % time range to display in feature view, in seconds
        umPerPix = 20;                    %

        % preview GUI params
        nLoadsMaxPreview = 30;            % number of time segments to load for preview
        nSecsLoadPreview = 1;           % recording duration per continuous segment to preview (in sec)

        % parameters for estimating firing rate
        frPeriod = 2;            % time period to determine the firing rate
        frFilterShape = 'triangle'; % {'triangle', 'rectangle'} kernel shape for temporal averaging
        frSampleRate = 1000;        % Resampled rate for the firing rate

        % aux-file parameters
        auxChan;                            % aux channel # to correlate with the unit firing rate
        auxFile = '';                       % aux channel file
        auxLabel = '';                      % label for the aux channel
        auxSampleRate = [];                 % sampling rate for aux file
        auxScale = 1;               		% scale factor for aux input
        nClustersShowAux = 10;              %

        % trial parameters
        trialFile = '';                     % .mat or .csv file containing timestamp in seconds unit. use any variable name.
        psthTimeLimits = [-1 5];                 % Time range to display PSTH (in seconds)
        psthTimeBin = .01;                    % Time bin for the PSTH histogram (in seconds)
        psthXTick = .2;                    % PSTH time tick mark spacing
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
        % autoMergeBy/autoMergeCriterion
        function set.autoMergeBy(obj, am)
            legalTypes = {'pearson', 'dist'};
            failMsg = sprintf('legal autoMergeBys are %s', strjoin(legalTypes, ', '));
            assert(sum(strcmp(am, legalTypes)) == 1, failMsg);
            obj.autoMergeBy = am;
        end
        function am = get.autoMergeCriterion(obj)
            obj.logOldP('autoMergeCriterion');
            am = obj.autoMergeBy;
        end
        function set.autoMergeCriterion(obj, am)
            obj.logOldP('autoMergeCriterion');
            obj.autoMergeBy = am;
        end

        % auxChan/iChan_aux
        function set.auxChan(obj, ac)
            assert(isempty(ac) || jrclust.utils.ismatrixnum(ac) && all(ac > 0), 'malformed auxChan');
            obj.auxChan = ac;
        end
        function ac = get.iChan_aux(obj)
            obj.logOldP('iChan_aux');
            ac = obj.auxChan;
        end
        function set.iChan_aux(obj, ac)
            obj.logOldP('iChan_aux');
            obj.auxChan = ac;
        end

        % auxFile/vcFile_aux
        function set.auxFile(obj, af)
            if isempty(af)
                obj.auxFile = '';
            else
                af_ = jrclust.utils.absPath(af);
                assert(isfile(af_), 'could not find aux file ''%s''', af);
                obj.auxFile = af_;
            end
        end
        function af = get.vcFile_aux(obj)
            obj.logOldP('vcFile_aux');
            af = obj.auxFile;
        end
        function set.vcFile_aux(obj, af)
            obj.logOldP('vcFile_aux');
            obj.auxFile = af;
        end

        % auxLabel/vcLabel_aux
        function set.auxLabel(obj, al)
            failMsg = 'auxLabel must be a string';
            assert(ischar(al), failMsg);
            obj.auxLabel = al;
        end
        function al = get.vcLabel_aux(obj)
            obj.logOldP('vcLabel_aux');
            al = obj.auxLabel;
        end
        function set.vcLabel_aux(obj, al)
            obj.logOldP('vcLabel_aux');
            obj.auxLabel = al;
        end

        % auxSampleRate/sRateHz_aux
        function set.auxSampleRate(obj, ar)
            failMsg = 'auxSampleRate must be a positive integer';
            assert(isempty(ar) || jrclust.utils.isscalarnum(ar) && ar == round(ar) && ar > 0, failMsg);
            obj.auxSampleRate = ar;
        end
        function ar = get.sRateHz_aux(obj)
            obj.logOldP('sRateHz_aux');
            ar = obj.auxSampleRate;
        end
        function set.sRateHz_aux(obj, ar)
            obj.logOldP('sRateHz_aux');
            obj.auxSampleRate = ar;
        end

        % auxScale/vrScale_aux
        function set.auxScale(obj, as)
            failMsg = 'auxScale must be a positive scalar';
            assert(jrclust.utils.isscalarnum(as) && as > 0, failMsg);
            assert(jrclust.utils.isscalarnum(as) && as > 0, failMsg);
            obj.auxScale = as;
        end
        function as = get.vrScale_aux(obj)
            obj.logOldP('vrScale_aux');
            as = obj.auxScale;
        end
        function set.vrScale_aux(obj, as)
            obj.logOldP('vrScale_aux');
            obj.auxScale = as;
        end

        % bitScaling/uV_per_bit
        function set.bitScaling(obj, bs)
            assert(jrclust.utils.isscalarnum(bs) && bs > 0, 'bad bitScaling factor');
            obj.bitScaling = bs;
        end
        function bs = get.uV_per_bit(obj)
            obj.logOldP('uV_per_bit');
            bs = obj.bitScaling;
        end
        function set.uV_per_bit(obj, bs)
            obj.logOldP('uV_per_bit');
            obj.bitScaling = bs;
        end

        % blankThresh/blank_thresh
        function set.blankThresh(obj, bt)
            assert((jrclust.utils.isscalarnum(bt) && bt >= 0) || (isnumeric(bt) && isempty(bt)), 'bad blankThresh');
            obj.blankThresh = bt;
        end
        function bt = get.blank_thresh(obj)
            obj.logOldP('blank_thresh');
            bt = obj.blankThresh;
        end
        function set.blank_thresh(obj, bt)
            obj.logOldP('blank_thresh');
            obj.blankThresh = bt;
        end

        % bytesPerSample
        function bp = get.bytesPerSample(obj)
            bp = jrclust.utils.typeBytes(obj.dataType);
        end

        % CARMode/vcCommonRef
        function set.CARMode(obj, cm)
            legalTypes = {'mean', 'median', 'none'};
            failMsg = sprintf('legal carModes are %s', strjoin(legalTypes, ', '));
            assert(ismember(cm, legalTypes), failMsg);
            obj.CARMode = cm;
        end
        function cm = get.vcCommonRef(obj)
            obj.logOldP('vcCommonRef');
            cm = obj.CARMode;
        end
        function set.vcCommonRef(obj, cm)
            obj.logOldP('vcCommonRef');
            obj.CARMode = cm;
        end

        % clusterFeature/vcFet
        function set.clusterFeature(obj, cf)
            if strcmp(cf, 'gpca')
                cf = 'pca';
                warning('gpca feature temporarily disabled; using pca');
            end
            legalTypes = {'cov', 'vpp', 'vmin', 'vminmax', 'energy', 'pca'};
            failMsg = sprintf('legal clusterFeatures are: %s', strjoin(legalTypes, ', '));
            assert(sum(strcmp(cf, legalTypes)) == 1, failMsg);
            obj.clusterFeature = cf;
        end
        function cf = get.vcFet(obj)
            obj.logOldP('vcFet');
            cf = obj.clusterFeature;
        end
        function set.vcFet(obj, cf)
            obj.logOldP('vcFet');
            obj.clusterFeature = cf;
        end

        % configFile/vcFile_prm
        function set.configFile(obj, cf)
            cf_ = jrclust.utils.absPath(cf);
            assert(isfile(cf_), 'could not find file ''%s''', cf);
            obj.configFile = cf_;
        end
        function cf = get.vcFile_prm(obj)
            obj.logOldP('vcFile_prm');
            cf = obj.configFile;
        end
        function set.vcFile_prm(obj, cf)
            obj.logOldP('vcFile_prm');
            obj.configFile = cf;
        end

        % dispFilter/vcFilter_show
        function set.dispFilter(obj, ft)
            legalTypes = {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'none'};
            failMsg = sprintf('legal filterTypes are: %s or ''''', strjoin(legalTypes, ', '));
            assert((ischar(ft) && isempty(ft)) || sum(strcmp(ft, legalTypes)) == 1, failMsg);
            obj.dispFilter = ft;
        end
        function ft = get.vcFilter_show(obj)
            obj.logOldP('vcFilter_show');
            ft = obj.dispFilter;
        end
        function set.vcFilter_show(obj, ft)
            obj.logOldP('vcFilter_show');
            obj.dispFilter = ft;
        end

        % dispFeature/vcFet_show
        function set.dispFeature(obj, df)
            % canonicalize synonyms
            if strcmp(df, 'vmin')
                df = 'vpp';
            elseif strcmp(df, 'gpca')
                df = 'pca';
            elseif strcmp(df, 'private pca')
                df = 'ppca';
            elseif strcmp(df, 'spacetime')
                df = 'cov';
            end

            legalTypes = {'cov', 'kilosort', 'pca', 'ppca', 'vpp'};
            failMsg = sprintf('legal dispFeatures are %s', strjoin(legalTypes, ', '));
            assert(ismember(df, legalTypes), failMsg);

            if ~strcmp(df, 'pca') % reset pcPair
                obj.pcPair = [1 2]; %#ok<MCSUP>
            end
            obj.dispFeature = df;
        end
        function df = get.vcFet_show(obj)
            obj.logOldP('vcFet_show');
            df = obj.dispFeature;
        end
        function set.vcFet_show(obj, df)
            obj.logOldP('vcFet_show');
            obj.dispFeature = df;
        end

        % dispTimeLimits/tlim
        function set.dispTimeLimits(obj, tl)
            assert((isempty(tl) && isnumeric(tl)) || (jrclust.utils.ismatrixnum(tl) && all(tl >= 0)), 'bad dispTimeLimits');
            if numel(tl) == 1
                tl = [0 tl];
            end
            assert(all(size(tl) == [1 2]) && tl(1) < tl(2), 'bad dispTimeLimits');
            obj.dispTimeLimits = tl;
        end
        function tl = get.tlim(obj)
            obj.logOldP('tlim');
            tl = obj.dispTimeLimits;
        end
        function set.tlim(obj, tl)
            obj.logOldP('tlim');
            obj.dispTimeLimits = tl;
        end

        % dataType/vcDataType
        function set.dataType(obj, dt)
            legalTypes = {'int16', 'uint16', 'int32', 'uint32', 'single', 'double'};
            failMsg = sprintf('legal dtypes are: %s', strjoin(legalTypes, ', '));
            assert(ismember(dt, legalTypes), failMsg);
            obj.dataType = dt;
        end
        function dt = get.vcDataType(obj)
            obj.logOldP('vcDataType');
            dt = obj.dataType;
        end
        function set.vcDataType(obj, dt)
            obj.logOldP('vcDataType');
            obj.dataType = dt;
        end

        % evtDetectRad/maxDist_site_spk_um
        function set.evtDetectRad(obj, ed)
            assert(jrclust.utils.isscalarnum(ed) && ed > 0, 'evtDetectRad must be a positive scalar');
            obj.evtDetectRad = ed;
        end
        function ed = get.maxDist_site_spk_um(obj)
            obj.logOldP('maxDist_site_spk_um');
            ed = obj.evtDetectRad;
        end
        function set.maxDist_site_spk_um(obj, ed)
            obj.logOldP('maxDist_site_spk_um');
            obj.evtDetectRad = ed;
        end

        % evtManualThreshSamp/spkThresh
        function mt = get.evtManualThreshSamp(obj)
            mt = obj.evtManualThresh / obj.bitScaling;
        end
        function mt = get.spkThresh(obj)
            obj.logOldP('spkThresh');
            mt = obj.evtManualThreshSamp;
        end

        % evtManualThresh/spkThresh_uV
        function set.evtManualThresh(obj, mt)
            assert(isempty(mt) || (jrclust.utils.isscalarnum(mt) && mt ~= 0)); % TODO: positive or negative?
            obj.evtManualThresh = mt;
        end
        function mt = get.spkThresh_uV(obj)
            obj.logOldP('spkThresh_uV');
            mt = obj.evtManualThresh;
        end
        function set.spkThresh_uV(obj, mt)
            obj.logOldP('spkThresh_uV');
            obj.evtManualThresh = mt;
        end

        % evtMergeRad/maxDist_site_um
        function set.evtMergeRad(obj, em)
            assert(jrclust.utils.isscalarnum(em) && em >= 0, 'evtMergeRad must be a nonnegative scalar');
            obj.evtMergeRad = em;
        end
        function em = get.maxDist_site_um(obj)
            obj.logOldP('maxDist_site_um');
            em = obj.evtMergeRad;
        end
        function set.maxDist_site_um(obj, em)
            obj.logOldP('maxDist_site_um');
            obj.evtMergeRad = em;
        end

        % evtWindow/spkLim_ms
        function set.evtWindow(obj, ew)
            assert(ismatrix(ew) && all(size(ew) == [1 2]) && ew(1) < 0 && ew(2) > 0, 'degenerate evtWindow');
            obj.evtWindow = ew;
        end
        function ew = get.spkLim_ms(obj)
            obj.logOldP('spkLim_ms');
            ew = obj.evtWindow;
        end
        function set.spkLim_ms(obj, ew)
            obj.logOldP('spkLim_ms');
            obj.evtWindow = ew;
        end

        % evtWindowRaw/spkLim_raw_ms
        function set.evtWindowRaw(obj, ew)
            assert(ismatrix(ew) && all(size(ew) == [1 2]) && ew(1) < 0 && ew(2) > 0, 'degenerate evtWindowRaw');
            obj.evtWindowRaw = ew;
        end
        function ew = get.spkLim_raw_ms(obj)
            obj.logOldP('spkLim_raw_ms');
            ew = obj.evtWindowRaw;
        end
        function set.spkLim_raw_ms(obj, ew)
            obj.logOldP('spkLim_raw_ms');
            obj.evtWindowRaw = ew;
        end

        % evtWindowRawSamp/spkLim_raw
        function ew = get.evtWindowRawSamp(obj)
            ew = round(obj.evtWindowRaw * obj.sampleRate / 1000);
        end
        function set.evtWindowRawSamp(obj, ew)
            obj.evtWindowRaw = ew * 1000 / obj.sampleRate;
        end
        function ew = get.spkLim_raw(obj)
            obj.logOldP('spkLim_raw');
            ew = obj.evtWindowRawSamp;
        end
        function set.spkLim_raw(obj, ew)
            obj.logOldP('spkLim_raw');
            obj.evtWindowRawSamp = ew;
        end

        % evtWindowSamp/spkLim
        function ew = get.evtWindowSamp(obj)
            ew = round(obj.evtWindow * obj.sampleRate / 1000);
        end
        function set.evtWindowSamp(obj, ew)
            obj.evtWindow = ew * 1000 / obj.sampleRate;
        end
        function ew = get.spkLim(obj)
            obj.logOldP('spkLim');
            ew = obj.evtWindowSamp;
        end
        function set.spkLim(obj, ew)
            obj.logOldP('spkLim');
            obj.evtWindowSamp = ew;
        end

        % fftThresh/fft_thresh
        function set.fftThresh(obj, ft)
            assert(jrclust.utils.isscalarnum(ft) && ft >= 0, 'fftThresh must be a nonnegative scalar');
            obj.fftThresh = ft;
        end
        function ft = get.fft_thresh(obj)
            obj.logOldP('fft_thresh');
            ft = obj.fftThresh;
        end
        function set.fft_thresh(obj, ft)
            obj.logOldP('fft_thresh');
            obj.fftThresh = ft;
        end

        % filtOrder
        function set.filtOrder(obj, fo)
            assert(jrclust.utils.isscalarnum(fo) && fo > 0, 'bad filtOrder');
            obj.filtOrder = fo;
        end

        % filterType/vcFilter
        function set.filterType(obj, ft)
            legalTypes = {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'none'};
            assert(sum(strcmp(ft, legalTypes)) == 1, 'legal filterTypes are: %s', strjoin(legalTypes, ', '));
            obj.filterType = ft;
        end
        function ft = get.vcFilter(obj)
            obj.logOldP('vcFilter');
            ft = obj.filterType;
        end
        function set.vcFilter(obj, ft)
            obj.logOldP('vcFilter');
            obj.filterType = ft;
        end

        % frFilterShape/filter_shape_rate
        function set.frFilterShape(obj, fr)
            legalTypes = {'triangle', 'rectangle'};
            failMsg = sprintf('legal frFilterShapes are %s', strjoin(legalTypes, ', '));
            assert(ismember(fr, legalTypes), failMsg);
            obj.frFilterShape = fr;
        end
        function fr = get.filter_shape_rate(obj)
            obj.logOldP('filter_shape_rate');
            fr = obj.frFilterShape;
        end
        function set.filter_shape_rate(obj, fr)
            obj.logOldP('filter_shape_rate');
            obj.frFilterShape = fr;
        end

        % frPeriod/filter_sec_rate
        function set.frPeriod(obj, fr)
            failMsg = 'frPeriod must be a positive scalar';
            assert(jrclust.utils.isscalarnum(fr) && fr > 0, failMsg);
            obj.frPeriod = fr;
        end
        function fr = get.filter_sec_rate(obj)
            obj.logOldP('filter_sec_rate');
            fr = obj.frPeriod;
        end
        function set.filter_sec_rate(obj, fr)
            obj.logOldP('filter_sec_rate');
            obj.frPeriod = fr;
        end

        % frSampleRate/sRateHz_rate
        function set.frSampleRate(obj, fr)
            failMsg = 'frSampleRate must be a positive integer';
            assert(jrclust.utils.isscalarnum(fr) && fr == round(fr) && fr > 0, failMsg);
            obj.frSampleRate = fr;
        end
        function fr = get.sRateHz_rate(obj)
            obj.logOldP('sRateHz_rate');
            fr = obj.frSampleRate;
        end
        function set.sRateHz_rate(obj, fr)
            obj.logOldP('sRateHz_rate');
            obj.frSampleRate = fr;
        end

        % freqLimBP
        function set.freqLimBP(obj, fl)
            assert(jrclust.utils.ismatrixnum(fl) && all(size(fl) == [1 2]) && all(fl >= 0), 'bad freqLimBP');
            obj.freqLimBP = fl;
        end

        % gainBoost/gain_boost
        function set.gainBoost(obj, gb)
            assert(jrclust.utils.isscalarnum(gb) && gb > 0, 'gainBoost must be a positive scalar');
            obj.gainBoost = gb;
        end
        function gb = get.gain_boost(obj)
            obj.logOldP('gain_boost');
            gb = obj.gainBoost;
        end
        function set.gain_boost(obj, gb)
            obj.logOldP('gain_boost');
            obj.gainBoost = gb;
        end

        % gtFile/vcFile_gt
        function set.gtFile(obj, gf)
            if isempty(gf)
                obj.gtFile = '';
            else
                gf_ = jrclust.utils.absPath(gf);
                assert(isfile(gf_), 'could not find ground-truth file ''%s''', gf);
                obj.gtFile = gf_;
            end
        end
        function gf = get.vcFile_gt(obj)
            obj.logOldP('vcFile_gt');
            gf = obj.gtFile;
        end
        function set.vcFile_gt(obj, gf)
            obj.logOldP('vcFile_gt');
            obj.gtFile = gf;
        end

        % headerOffset/header_offset
        function set.headerOffset(obj, ho)
            assert(jrclust.utils.isscalarnum(ho) && ho >= 0, 'invalid headerOffset');
            obj.headerOffset = ho;
        end
        function ho = get.header_offset(obj)
            obj.logOldP('header_offset');
            ho = obj.headerOffset;
        end
        function set.header_offset(obj, ho)
            obj.logOldP('header_offset');
            obj.headerOffset = ho;
        end

        % ignoreSites/viSiteZero
        function set.ignoreSites(obj, ig)
            assert(jrclust.utils.ismatrixnum(ig) && all(ig > 0), 'degenerate ignoreSites');
            % don't manually ignore sites that are automatically ignored
            obj.ignoreSites = ig;
        end
        function ig = get.viSiteZero(obj)
            obj.logOldP('viSiteZero');
            ig = obj.ignoreSites;
        end
        function set.viSiteZero(obj, ig)
            obj.logOldP('viSiteZero');
            obj.ignoreSites = ig;
        end

        % lfpSampleRate/sRateHz_lfp
        function set.lfpSampleRate(obj, lf)
            assert(jrclust.utils.isscalarnum(lf) && lf > 0, 'bad lfpSampleRate');
            obj.lfpSampleRate = lf;
        end
        function lf = get.sRateHz_lfp(obj)
            obj.logOldP('sRateHz_lfp');
            lf = obj.lfpSampleRate;
        end
        function set.sRateHz_lfp(obj, lf)
            obj.logOldP('sRateHz_lfp');
            obj.lfpSampleRate = lf;
        end

        % loadTimeLimits/tlim_load
        function set.loadTimeLimits(obj, tl)
            assert((isempty(tl) && isnumeric(tl)) || ...
                   (jrclust.utils.ismatrixnum(tl) && all(size(tl) == [1 2]) && all(tl > 0) && tl(1) < tl(2)), 'bad loadTimeLimits');
            obj.loadTimeLimits = tl;
        end
        function tl = get.tlim_load(obj)
            obj.logOldP('tlim_load');
            tl = obj.loadTimeLimits;
        end
        function set.tlim_load(obj, tl)
            obj.logOldP('tlim_load');
            obj.loadTimeLimits = tl;
        end

        % log10DeltaCut/delta1_cut
        function set.log10DeltaCut(obj, dc)
            assert(jrclust.utils.isscalarnum(dc), 'log10DeltaCut must be a numeric scalar');
            obj.log10DeltaCut = dc;
        end
        function dc = get.delta1_cut(obj)
            obj.logOldP('delta1_cut');
            dc = obj.log10DeltaCut;
        end
        function set.delta1_cut(obj, dc)
            obj.logOldP('delta1_cut');
            obj.log10DeltaCut = dc;
        end

        % log10RhoCut/rho_cut
        function set.log10RhoCut(obj, rc)
            assert(jrclust.utils.isscalarnum(rc), 'log10RhoCut must be a numeric scalar');
            obj.log10RhoCut = rc;
        end
        function rc = get.rho_cut(obj)
            obj.logOldP('rho_cut');
            rc = obj.log10RhoCut;
        end
        function set.rho_cut(obj, rc)
            obj.logOldP('rho_cut');
            obj.log10RhoCut = rc;
        end

        % maxBytesLoad/MAX_BYTES_LOAD
        function set.maxBytesLoad(obj, mb)
            assert(isempty(mb) || jrclust.utils.isscalarnum(mb) && mb > 0, 'maxBytesLoad must be a positive scalar');
            obj.maxBytesLoad = mb;
        end
        function mb = get.MAX_BYTES_LOAD(obj)
            obj.logOldP('MAX_BYTES_LOAD');
            mb = obj.maxBytesLoad;
        end
        function set.MAX_BYTES_LOAD(obj, mb)
            obj.logOldP('MAX_BYTES_LOAD');
            obj.maxBytesLoad = mb;
        end

        % maxClustersSite/maxCluPerSite
        function set.maxClustersSite(obj, mc)
            failMsg = 'maxClustersSite must be a nonnegative integer';
            assert(jrclust.utils.isscalarnum(mc) && round(mc) == mc && mc > 0, failMsg);
            obj.maxClustersSite = mc;
        end
        function mc = get.maxCluPerSite(obj)
            obj.logOldP('maxCluPerSite');
            mc = obj.maxClustersSite;
        end
        function set.maxCluPerSite(obj, mc)
            obj.logOldP('maxCluPerSite');
            obj.maxClustersSite = mc;
        end

        % maxSecLoad/MAX_LOAD_SEC
        function set.maxSecLoad(obj, ms)
            assert(isempty(ms) || (jrclust.utils.isscalarnum(ms) && ms > 0), 'maxSecLoad must be a positive scalar');
            obj.maxSecLoad = ms;
        end
        function ms = get.MAX_LOAD_SEC(obj)
            obj.logOldP('MAX_LOAD_SEC');
            ms = obj.maxSecLoad;
        end
        function set.MAX_LOAD_SEC(obj, ms)
            obj.logOldP('MAX_LOAD_SEC');
            obj.maxSecLoad = ms;
        end

        % maxUnitSim
        function set.maxUnitSim(obj, mw)
            failMsg = 'maxUnitSim must be between 0 and 1';
            assert(jrclust.utils.isscalarnum(mw) && mw >= 0 && mw <= 1, failMsg);
            obj.maxUnitSim = mw;
        end

        % minClusterSize/min_count
        function set.minClusterSize(obj, mc)
            assert(jrclust.utils.isscalarnum(mc) && mc == round(mc) && mc > 0, 'minClusterSize must be a positive integer-valued scalar');
            obj.minClusterSize = mc;
        end
        function mc = get.min_count(obj)
            obj.logOldP('min_count');
            mc = obj.minClusterSize;
        end
        function set.min_count(obj, mc)
            obj.logOldP('min_count');
            obj.minClusterSize = mc;
        end

        % csFile_merge
        function mr = get.csFile_merge(obj)
            obj.logOldP('csFile_merge');
            mr = obj.rawRecordings;
        end

        % nClustersShowAux/nClu_show_aux
        function set.nClustersShowAux(obj, nc)
            failMsg = 'nClustersShowAux must be a positive integer';
            assert(jrclust.utils.isscalarnum(nc) && nc == round(nc) && nc > 0, failMsg);
            obj.nClustersShowAux = nc;
        end
        function nc = get.nClu_show_aux(obj)
            obj.logOldP('nClu_show_aux');
            nc = obj.nClustersShowAux;
        end
        function set.nClu_show_aux(obj, nc)
            obj.logOldP('nClu_show_aux');
            obj.nClustersShowAux = nc;
        end

        % nPeaksFeatures
        function set.nPeaksFeatures(obj, nf)
            assert(jrclust.utils.isscalarnum(nf) && ~isempty(intersect(nf, [1 2 3])), 'nPeaksFeatures must be 1, 2, or 3');
            obj.nPeaksFeatures = nf;
        end

        % nPassesMerge/nRepeat_merge
        function set.nPassesMerge(obj, np)
            failMsg = 'nPassesMerge must be a nonnegative integer';
            assert(jrclust.utils.isscalarnum(np) && np == round(np) && np >= 0, failMsg);
            obj.nPassesMerge = np;
        end
        function np = get.nRepeat_merge(obj)
            obj.logOldP('nRepeat_merge');
            np = obj.nPassesMerge;
        end
        function set.nRepeat_merge(obj, np)
            obj.logOldP('nRepeat_merge');
            obj.nPassesMerge = np;
        end

        % nSamplesPad/nPad_filt
        function set.nSamplesPad(obj, ns)
            assert(jrclust.utils.isscalarnum(ns) && ns >= 0, 'nSamplesPad must be a nonnegative scalar');
            obj.nSamplesPad = ns;
        end
        function ns = get.nPad_filt(obj)
            obj.logOldP('nPad_filt');
            ns = obj.nSamplesPad;
        end
        function set.nPad_filt(obj, ns)
            obj.logOldP('nPad_filt');
            obj.nSamplesPad = ns;
        end

        % nSiteDir/maxSite
        function set.nSiteDir(obj, ns)
            assert(jrclust.utils.ismatrixnum(ns) && (isempty(ns) || (jrclust.utils.isscalarnum(ns) && ns >= 0)), 'invalid value for nSiteDir');
            obj.nSiteDir = ns;
        end
        function ns = get.maxSite(obj)
            obj.logOldP('maxSite');
            ns = obj.nSiteDir;
        end
        function set.maxSite(obj, ns)
            obj.logOldP('maxSite');
            obj.nSiteDir = ns;
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
        function set.nSitesExcl(obj, ns)
            assert(jrclust.utils.ismatrixnum(ns) && (isempty(ns) || (jrclust.utils.isscalarnum(ns) && ns >= 0)), 'invalid value for nSitesExcl');
            obj.nSitesExcl = ns;
        end
        function ns = get.nSites_ref(obj)
            obj.logOldP('nSites_ref');
            ns = obj.nSitesExcl;
        end
        function set.nSites_ref(obj, ns)
            obj.logOldP('nSites_ref');
            obj.nSitesExcl = ns;
        end

        % outlierThresh/thresh_mad_clu
        function set.outlierThresh(obj, ot)
            failMsg = 'outlierThresh must be a nonnegative number';
            assert(jrclust.utils.isscalarnum(ot) && ot >= 0, failMsg);
            obj.outlierThresh = ot;
        end
        function ot = get.thresh_mad_clu(obj)
            obj.logOldP('thresh_mad_clu');
            ot = obj.outlierThresh;
        end
        function set.thresh_mad_clu(obj, ot)
            obj.logOldP('thresh_mad_clu');
            obj.outlierThresh = ot;
        end

        % outputDir
        function od = get.outputDir(obj)
            if isempty(obj.outputDir)
                od = fileparts(obj.configFile);
            else
                od = obj.outputDir;
            end
        end
        function set.outputDir(obj, od)
            failMsg = 'outputDir must be a directory';
            od_ = jrclust.utils.absPath(od);
            assert(isempty(od) || isdir(od_), failMsg);
            obj.outputDir = od_;
        end

        % probeFile/probe_file
        function set.probeFile(obj, pf)
            % first check where the config file is located
            pf_ = jrclust.utils.absPath(pf);
            % if we can't find it there, try the standard location
            if isempty(pf_)
                pf_ = jrclust.utils.absPath(pf, fullfile(jrclust.utils.basedir(), 'probes'));
            end
            assert(isfile(pf_), 'could not find probe file ''%s''', pf);
            obj.probeFile = pf_;
        end
        function pf = get.probe_file(obj)
            obj.logOldP('probe_file');
            pf = obj.probeFile;
        end
        function set.probe_file(obj, pf)
            obj.logOldP('probe_file');
            obj.probeFile = pf;
        end

        % probePad/vrSiteHW
        function set.probePad(obj, pp)
            assert(jrclust.utils.ismatrixnum(pp) && all(pp(:) > 0), 'malformed probePad');
            obj.probePad = pp;
        end
        function pp = get.vrSiteHW(obj)
            obj.logOldP('vrSiteHW');
            pp = obj.probePad;
        end
        function set.vrSiteHW(obj, pp)
            obj.logOldP('vrSiteHW');
            obj.probePad = pp;
        end

        % ramToGPUFactor/nLoads_gpu
        function set.ramToGPUFactor(obj, rg)
            assert(jrclust.utils.isscalarnum(rg) && rg > 0, 'ramToGPUFactor must be a positive scalar');
            obj.ramToGPUFactor = rg;
        end
        function rg = get.nLoads_gpu(obj)
            obj.logOldP('nLoads_gpu');
            rg = obj.ramToGPUFactor;
        end
        function set.nLoads_gpu(obj, rg)
            obj.logOldP('nLoads_gpu');
            obj.ramToGPUFactor = rg;
        end

        % randomSeed
        function set.randomSeed(obj, rs)
            failMsg = 'randomSeed must be a nonnegative integer';
            assert(jrclust.utils.isscalarnum(rs) && rs == round(rs) && rs >= 0, failMsg);
            obj.randomSeed = rs;
        end

        % refracInt/spkRefrac_ms
        function set.refracInt(obj, ri)
            assert(jrclust.utils.isscalarnum(ri) && ri > 0, 'refractory interval must be a positive scalar');
            obj.refracInt = ri;
        end
        function ri = get.spkRefrac_ms(obj)
            obj.logOldP('spkRefrac_ms');
            ri = obj.refracInt;
        end
        function set.spkRefrac_ms(obj, ri)
            obj.logOldP('spkRefrac_ms');
            obj.refracInt = ri;
        end

        % refracIntSamp/spkRefrac
        function ri = get.refracIntSamp(obj)
            ri = round(obj.refracInt * obj.sampleRate / 1000);
        end
        function set.refracIntSamp(obj, ri)
            obj.refracInt = ri * 1000 / obj.sampleRate;
        end
        function ri = get.spkRefrac(obj)
            obj.logOldP('spkRefrac');
            ri = obj.refracIntSamp;
        end
        function set.spkRefrac(obj, ri)
            obj.logOldP('spkRefrac');
            obj.refracIntSamp = ri;
        end

        % RDDetrendMode/vcDetrend_postclu
        function set.RDDetrendMode(obj, dm)
            legalTypes = {'global', 'local', 'logz', 'hidehiko', 'none'};
            failMsg = sprintf('legal RDDetrendModes are %s', strjoin(legalTypes, ', '));
            assert(sum(strcmp(dm, legalTypes)) == 1, failMsg);
            obj.RDDetrendMode = dm;
        end
        function dm = get.vcDetrend_postclu(obj)
            obj.logOldP('vcDetrend_postclu');
            dm = obj.RDDetrendMode;
        end
        function set.vcDetrend_postclu(obj, dm)
            obj.logOldP('vcDetrend_postclu');
            obj.RDDetrendMode = dm;
        end

        % sampleRate/sRateHz
        function set.sampleRate(obj, sr)
            assert(jrclust.utils.isscalarnum(sr) && sr > 0, 'bad sample rate');
            obj.sampleRate = sr;
        end
        function sr = get.sRateHz(obj)
            obj.logOldP('sRateHz');
            sr = obj.sampleRate;
        end
        function set.sRateHz(obj, sr)
            obj.logOldP('sRateHz');
            obj.sampleRate = sr;
        end

        % sessionName
        function sn = get.sessionName(obj)
            [~, sn, ~] = fileparts(obj.configFile);
        end

        % shankMap/viShank_site
        function set.shankMap(obj, sm)
            assert(jrclust.utils.ismatrixnum(sm) && all(sm > 0), 'malformed shankMap');
            obj.shankMap = sm;
        end
        function sm = get.viShank_site(obj)
            obj.logOldP('viShank_site');
            sm = obj.shankMap;
        end
        function set.viShank_site(obj, sm)
            obj.logOldP('viShank_site');
            obj.shankMap = sm;
        end

        % showRaw/fWav_raw_show
        function set.showRaw(obj, sr)
            obj.showRaw = 1 && sr;
        end
        function sr = get.fWav_raw_show(obj)
            obj.logOldP('fWav_raw_show');
            sr = obj.showRaw;
        end
        function set.fWav_raw_show(obj, sr)
            obj.logOldP('fWav_raw_show');
            obj.showRaw = sr;
        end

        % siteCorrThresh/thresh_corr_bad_site
        function set.siteCorrThresh(obj, ct)
            assert(jrclust.utils.isscalarnum(ct) && ct >= 0 && ct < 1, 'bad siteCorrThresh');
            obj.siteCorrThresh = ct;
        end
        function ct = get.thresh_corr_bad_site(obj)
            obj.logOldP('thresh_corr_bad_site');
            ct = obj.siteCorrThresh;
        end
        function set.thresh_corr_bad_site(obj, ct)
            obj.logOldP('thresh_corr_bad_site');
            obj.siteCorrThresh = ct;
        end

        % siteLoc/mrSiteXY
        function set.siteLoc(obj, pg)
            assert(jrclust.utils.ismatrixnum(pg) && all(pg(:) >= 0), 'malformed siteLoc');
            obj.siteLoc = pg;
        end
        function pg = get.mrSiteXY(obj)
            obj.logOldP('mrSiteXY');
            pg = obj.siteLoc;
        end
        function set.mrSiteXY(obj, pg)
            obj.logOldP('mrSiteXY');
            obj.siteLoc = pg;
        end

        % siteMap/viSite2Chan
        function set.siteMap(obj, cm)
            assert(jrclust.utils.ismatrixnum(cm) && all(cm > 0), 'malformed siteMap');
            obj.siteMap = cm;
        end
        function cm = get.viSite2Chan(obj)
            obj.logOldP('viSite2Chan');
            cm = obj.siteMap;
        end
        function set.viSite2Chan(obj, cm)
            obj.logOldP('viSite2Chan');
            obj.siteMap = cm;
        end

        % siteNeighbors/miSites
        function set.siteNeighbors(obj, sn) % danger zone: don't set this manually
            assert(jrclust.utils.ismatrixnum(sn), 'siteNeighbors must be a numeric matrix');
            obj.siteNeighbors = sn;
        end
        function sn = get.miSites(obj)
            obj.logOldP('miSites');
            sn = obj.siteNeighbors;
        end
        function set.miSites(obj, sn)
            obj.logOldP('miSites');
            obj.siteNeighbors = sn;
        end

        % threshFile/vcFile_thresh
        function set.threshFile(obj, tf)
            if isempty(tf)
                obj.threshFile = '';
            else
                failMsg = sprintf('could not find threshold file ''%s''', tf);
                tf_ = jrclust.utils.absPath(tf);
                assert(isfile(tf_), failMsg);
                obj.threshFile = tf_;
            end
        end
        function gf = get.vcFile_thresh(obj)
            obj.logOldP('vcFile_thresh');
            gf = obj.threshFile;
        end
        function set.vcFile_thresh(obj, gf)
            obj.logOldP('vcFile_thresh');
            obj.threshFile = gf;
        end

        % trialFile/vcFile_trial
        function set.trialFile(obj, tf)
            if isempty(tf)
                obj.trialFile = tf;
            else
                failMsg = sprintf('could not find trial file ''%s''', tf);
                tf_ = jrclust.utils.absPath(tf);
                assert(isfile(tf_), failMsg);
                obj.trialFile = tf_;
            end
        end
        function tf = get.vcFile_trial(obj)
            obj.logOldP('vcFile_trial');
            tf = obj.trialFile;
        end
        function set.vcFile_trial(obj, tf)
            obj.logOldP('vcFile_trial');
            obj.trialFile = tf;
        end

        % useGPU/fGpu
        function set.useGPU(obj, ug)
            % use GPU if and only if present, accessible, and asked for
            ug = ug && license('test', 'Distrib_Computing_Toolbox') && gpuDeviceCount() > 0;
            obj.useGPU = ug;
        end
        function ug = get.fGpu(obj)
            obj.logOldP('fGpu');
            ug = obj.useGPU;
        end
        function set.fGpu(obj, ug)
            obj.logOldP('fGpu');
            obj.useGPU = ug;
        end

        % userFiltKernel/vnFilter_user
        function set.userFiltKernel(obj, uf)
            assert(jrclust.utils.ismatrixnum(uf));
            obj.userFiltKernel = uf;
        end
        function uf = get.vnFilter_user(obj)
            obj.logOldP('vnFilter_user');
            uf = obj.userFiltKernel;
        end
        function set.vnFilter_user(obj, uf)
            obj.logOldP('vnFilter_user');
            obj.userFiltKernel = uf;
        end

        % verbose/fVerbose
        function vb = get.verbose(obj)
            vb = obj.verbose && ~obj.batchMode;
        end
        function set.verbose(obj, vb)
            obj.verbose = 1 && vb;
        end
        function vb = get.fVerbose(obj)
            obj.logOldP('fVerbose');
            vb = obj.verbose;
        end
        function set.fVerbose(obj, vb)
            obj.logOldP('fVerbose');
            obj.verbose = vb;
        end
    end
end
