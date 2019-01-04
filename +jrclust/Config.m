classdef Config < dynamicprops
    %CONFIG JRCLUST session configuration
    % replacement for P struct

    %% OBJECT-LEVEL PROPERTIES
    properties (Hidden, SetObservable, SetAccess=private)
        isLoaded;
        errMsg;
        isError;
        oldPcount;
    end

    %% OLD-STYLE PARAMS, publicly settable (will be deprecated after a grace period)
    properties (SetObservable, Dependent, Hidden, Transient)
        % cvrDepth_drift;           % doesn't appear to be used (was {})
        % maxSite_detect;           % doesn't appear to be used (was 1.5)
        % maxSite_dip;              % doesn't appear to be used (was [])
        % maxSite_fet;              % doesn't appear to be used (was [])
        % maxSite_merge;            % doesn't appear to be used (was [])
        % maxSite_pix;              % doesn't appear to be used (was 4.5)
        % maxSite_show;             % appears to be synonymous with nSiteDir/maxSite
        % maxSite_sort;             % doesn't appear to be used (was [])
        % rejectSpk_mean_thresh;    % appears to be synonymous with blankThresh/blank_thresh
        autoMergeCriterion;         % => autoMergeBy
        blank_thresh;               % => blankThresh
        csFile_merge;               % => multiRaw
        delta1_cut;                 % => log10DeltaCut
        fEllip;                     % => useElliptic
        fft_thresh;                 % => fftThreshMad
        fGpu;                       % => useGPU
        fImportKilosort;            % => fImportKsort
        fRepeat_clu;                % => repeatLower
        fVerbose;                   % => verbose
        fWav_raw_show;              % => showRaw
        gain_boost;                 % => gainBoost
        header_offset;              % => headerOffset
        MAX_BYTES_LOAD;             % => maxBytesLoad
        MAX_LOAD_SEC;               % => maxSecLoad
        maxCluPerSite;              % => maxClustersSite
        maxDist_site_um             % => evtMergeRad
        maxDist_site_spk_um;        % => evtDetectRad
        maxSite;                    % => nSiteDir
        min_count;                  % => minClusterSize
        mrSiteXY;                   % => siteLoc
        nLoads_gpu;                 % => ramToGPUFactor
        nPad_filt;                  % => nSamplesPad
        nRepeat_merge;              % => nPassesMerge
        nSites_ref;                 % => nSitesExcl
        nSkip_lfp;                  % => lfpDsFactor
        probe_file;                 % => probeFile
        rho_cut;                    % => log10RhoCut
        spkLim;                     % => evtWindowSamp
        spkLim_ms;                  % => evtWindowms
        spkLim_raw;                 % => evtWindowRawSamp
        spkLim_raw_factor;          % => evtWindowRawFactor
        spkLim_raw_ms;              % => evtWindowRawms
        spkRefrac;                  % => refracIntSamp
        spkRefrac_ms;               % => refracIntms
        spkThresh;                  % => evtManualThresh
        spkThresh_uV;               % => evtManualThreshuV
        sRateHz;                    % => sampleRate
        sRateHz_lfp;                % => lfpSampleRate
        thresh_corr_bad_site;       % => siteCorrThresh
        thresh_mad_clu;             % => outlierThresh
        tlim;                       % => dispTimeLimits
        tlim_load;                  % => loadTimeLimits
        uV_per_bit;                 % => bitScaling
        vcCommonRef;                % => carMode
        vcDataType;                 % => dtype
        vcDetrend_postclu;          % => rlDetrendMode
        vcFet;                      % => clusterFeature
        vcFet_show;                 % => dispFeature
        vcFile;                     % => singleRaw
        vcFile_gt;                  % => gtFile
        vcFile_prm;                 % => configFile
        vcFile_thresh;              % => threshFile
        vcFilter;                   % => filterType
        vcFilter_show;              % => dispFilter
        viChan_aux;                 % => auxSites
        viShank_site;               % => shankMap
        vnFilter_user;              % => userFiltKernel
        vrSiteHW;                   % => probePad
        viSite2Chan;                % => siteMap
        viSiteZero;                 % => ignoreSites
    end

    %% OLD-STLYE PARAMS, not publicly settable
    properties (Dependent, SetAccess=private, Hidden)
        miSites;                    % => siteNeighbors
    end
    
    %% NEW-STYLE PARAMS, publicly settable
    properties (SetObservable)
        % computation params
        gpuLoadFactor = 5;          % GPU memory usage factor (4x means 1/4 of GPU memory can be loaded)
        randomSeed = 0;             % random seed
        ramToGPUFactor = 8;         % ratio: RAM / (GPU memory) (increase this number if GPU memory error)
        useGPU = true;              % use GPU in computation if true
        verbose = true;             % be chatty while processing

        % file location params
        outputDir = '';             % directory in which to place output files

        % recording params
        auxSites;                   %
        bitScaling = 0.30518;       % bit scaling factor (uV/bit)
        configFile;                 % parameter file
        dtype = 'int16';            % raw data binary format
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
        useElliptic = true;         % use elliptic filter if true (and only if filterType='bandpass')
        fftThreshMAD = 0;           % automatically remove frequency outliers (unit:MAD, 10 recommended, 0 to disable). Verify by running "jrc traces" and press "p" to view the power spectrum.
        filtOrder = 3;              % bandpass filter order
        filterType = 'ndiff';       % filter to use {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'fftdiff', 'none'}
        freqLim = [300 3000];       % frequency cut-off limit for filterType='bandpass' (ignored otherwise)
        freqLimNotch = [];
        freqLimStop = [];
        loadTimeLimits = [];        % time range of recording to load, in s (use whole range if empty)
        maxBytesLoad = [];          % default memory loading block size (bytes)
        maxSecLoad = [];            % maximum loading duration (seconds) (overrides 'maxBytesLoad')
        ndist_filt = 5;             % undocumented
        nSamplesPad = 100;          % number of samples to overlap between multiple loading (filter edge safe)
        userFiltKernel = [];        % custom filter kernel (optional unless filterType='user')
        carMode = 'mean';           % common average referencing mode (one of 'none', 'mean', 'median', or 'whiten')

        % spike detection params
        blankThresh = [];           % reject spikes exceeding the channel mean after filtering (MAD unit), ignored if [] or 0
        evtDetectRad = 75;          % radius for extracting waveforms, in microns (used if nSiteDir and nSitesExcl are empty)
        evtManualThreshuV = [];     % manual spike detection threshold, in microvolts
        evtMergeRad = 50;           % radius of spike event merging, in microns
        evtWindowms = [-0.25 0.75]; % interval around event to extract filtered spike waveforms, in ms
        evtWindowRawms;             % interval around event to extract raw spike waveforms, in ms
        evtWindowRawFactor = 2;     % ratio of raw samples to filtered samples to extract if evtWindowRawms is not set
        fGroup_shank = false;       % group all sites in the same shank if true
        ignoreSites = [];           % sites to manually ignore in the sorting
        nDiff_filt = 2;             % Differentiation filter for filterType='sgdiff', ignored otherwise. Set to [] to disable. 2n+1 samples used for centered differentiation
        nneigh_min_detect = 0;      % Min. number of neighbors near the spike below threshold. choose between [0,1,2]
        nSiteDir;                   % number of neighboring sites to group in each direction (TODO: deprecate this)
        nSitesExcl;                 % number of sites to exclude from the spike waveform group
        refracIntms = 0.25;         % spike refractory interval, in ms
        siteCorrThresh = 0;         % reject bad sites based on max correlation with neighboring sites, using raw waveforms; ignored if 0
        spkThresh_max_uV = [];      % maximum absolute amp. allowed
        threshFile = '';            % path to .mat file storing spike detection thresholds (created by 'preview' GUI)

        % feature extraction params
        clusterFeature = 'pca';     % feature to use in clustering
        fInterp_fet = true;         % interpolate waveforms for feature projection to find optimal delay (2x interp) if true
        fSpatialMask_clu = false;   % apply spatial mask calculated from the distances between sites to the peak site (half-scale: evtDetectRad)
        min_sites_mask = 5;         % minimum number of sites to have to apply spatial mask
        nFet_use = 2;               % undocumented
        time_feature_factor;        % undocumented

        % clustering params
        autoMergeBy = 'xcorr';      % metric to use when automerging clusters
        dc_percent = 2;             % percentile at which to cut off distance in rho computation
        fDrift_merge = true;        % compute multiple waveforms at three drift locations based on the spike position if true
        log10DeltaCut = 0.6;        % the base-10 log of the delta cutoff value
        log10RhoCut = -2.5;         % the base-10 log of the rho cutoff value
        maxClustersSite = 20;       % maximum number of clusters per site if local detrending is used
        minClusterSize = 30;        % minimum cluster size (set to 2*#features if lower)
        nInterp_merge = 1;          % Interpolation factor for the mean unit waveforms, set to 1 to disable
        nPassesMerge = 10;          % number of passes for unit mean raw waveform-based merging
        outlierThresh = 7.5;        % threshold to remove outlier spikes for each cluster, in MAD
        nTime_clu = 1;              % number of time periods over which to cluster separately (later to be merged after clustering)
        repeatLower = false;        % repeat clustering for the bottom half of the cluster amplitudes if true
        rlDetrendMode = 'global';   % 
        spkLim_factor_merge = 1;    % Waveform range for computing the correlation. spkLim_factor_merge <= spkLim_raw_factor_merge. circa v3.1.8

        % display params
        dispFeature = 'vpp';        % feature to display in time/projection views
        dispFilter = '';            % 
        dispTimeLimits = [0 0.2];   % time range to display (in seconds)
        fText = true;               % 
        nShow = 200;                % maximum number of traces to show [D?# spikes to show]
        nShow_proj = 500;           % maximum number of features to show in projection
        nSitesFigProj = 5;          % number of sites to display in the feature projection view
        nTime_traces = 1;           % number of time segments to display. Set to 1 to show one continuous time segment
        nSpk_show = 30;             % show spike waveforms for manual clustering
        pcPair = [1 2];             % PC projection to show (1 vs 2; 1 vs 3; 2 vs 3), can be toggled
        showRaw = false;            % show raw waveforms in main view if true
        time_tick_show = [];        % 
        tLimFigProj = [];           % time range to display in feature view, in seconds
        um_per_pix = 20;            % 

        % to get to, eventually
        LineStyle = '';
        MAX_LOG = 5;
        S_imec3 = [];
        blank_period_ms = 5;
        corrLim = [0.9 1];
        cviShank = [];
        dc_factor = 1;
        dc_frac = [];
        dinput_imec_trial = 1;
        duration_file = [];
        fAddCommonRef = false;
        fAverageTrial_psth = true;
        fCacheRam = true;
        fCheckSites = false;
        fDetectBipolar = false;
        fDiscard_count = true;
        fInverse_file = false;
        fImportKsort = false;
        fLoad_lfp = false;
        fMeanSite = true;
        fMeanSiteRef = false;
        fMeanSite_drift = false;
        fMinNorm_wav = false;
        fNormRhoDelta = true;
        fParfor = true;
        fPcaDetect = false;
        fProcessEven = false;
        fProcessOdd = false;
        fProcessReverseOrder = false;
        fProj_sort = false;
        fRamCache = true;
        fRejectSpk_vpp = false;
        fRms_detect = false;
        fRun = true;
        fSaveEvt = true;
        fSavePlot_RD = true;
        fSaveRawSpk = false;
        fSaveSpk = true;
        fShowAllSites = false;
        fSingleColumn_track = true;
        fSmooth_track = true;
        fSpike_show = true;
        fTranspose_bin = true;
        fUseCache_track = false;
        fUseLfp_track = true;
        fWhiten_traces = false;
        filter_sec_rate = 2;
        filter_shape_rate = 'triangle';
        flim_vid = [];
        freqLimNotch_lfp = [];
        freqLim_corr = [15 150];
        freqLim_excl_track = [58 62];
        freqLim_lfp = [];
        freqLim_track = [15 150];
        iChan_aux = [];
        iChan_vid = [];
        iClu_show = [];
        iGpu = 1;
        load_fraction_track = [];
        maxAmp = 250;
        maxAmp_lfp = 1000;
        maxDist_site_merge_um = 35;
        maxLfpSdZ = 4.5;
        maxSite_track = [2 3 4 5 6 7 8];
        maxWavCor = 0.98;
        max_shift_track = [];
        mrColor_proj = [213 219 235; 0 130 196; 240 119 22]/256;
        nBytes_file = [];
        nClu_show_aux = 10;
        nLoads_max_preview = 30;
        nMinAmp_ms = 0;
        nPcPerChan = 1;
        nPc_dip = 3;
        nSites_excl_ref = 6;
        nSkip_show = 1;
        nSkip_whiten = 10;
        nSmooth_ms_psth = 50;
        nT_drift = [];
        nThreads = 128;
        nw_lcm_track = 1;
        offset_sec_preview = 0;
        pix_per_sec_track = [];
        qqFactor = 5;
        qqSample = 4;
        rateLim_psth = [];
        refrac_factor = 2;
        rms_filt_ms = 0;
        sRateHz_aux = [];
        sRateHz_rate = 1000;
        sec_per_load_preview = 1;
        slopeLim_ms = [0.05 0.35];
        spkLim_ms_fet = [-0.25 0.75];
        tBin_track = 9;
        tRefrac_trial = 0.001;
        tbin_drift = [];
        tbin_psth = 0.01;
        template_file = '';
        thresh_automerge_pca = [];
        thresh_corr_track = [];
        thresh_merge_clu = 0;
        thresh_sd_ref = 5;
        thresh_split_clu = 0;
        thresh_trial = [];
        tlim_clu = [];
        tlim_lfp = [0 5];
        tlim_psth = [-1 5];
        tlim_vid = [];
        vcCluDist = 'eucldist';
        vcCluWavMode = 'mean';
        vcDate_file = '';
        vcDc_clu = 'distr';
        vcFile_aux = '';
        vcFile_bonsai = '';
        vcFile_lfp = '';
        vcFile_trial = '';
        vcFile_vid = '';
        vcFilter_detect = '';
        vcLabel_aux = '';
        vcMode_track = 'mt_cpsd2_mr';
        vcSpatialFilter = 'none';
        vcSpkRef = 'nmean';
        viChan_bin = [];
        viChan_show = [];
        viDepth_excl_track = [];
        viDepth_track = [];
        viSite_bad_track = [];
        vrScale_aux = 1;
        xtick_psth = 0.2;
        ybin_drift = 2;
    end

    %% NEW-STYLE PARAMS, not publicly settable
    properties (SetAccess=private)
        siteNeighbors;              % indices of neighbors for each site
    end

    %% COMPUTED PARAMS
    properties (SetObservable, Dependent)
        bytesPerSample;             % byte count for each discrete sample
        evtManualThresh;            % evtManualThreshuV / bitScaling
        evtWindowRawSamp;           % interval around event to extract raw spike waveforms, in samples
        evtWindowSamp;              % interval around event to extract filtered spike waveforms, in samples
        nSitesEvt;                  % 2*nSiteDir + 1 - nSitesExcl
        refracIntSamp;              % spike refractory interval, in samples
        sessionName;                % name of prm file, without path or extensions
    end

    % here to ease the transition (use rawRecordings to access both)
    properties (Access=private, Hidden, SetObservable)
        multiRaw;                   % list of recording files to merge (empty if single file is used)
        singleRaw;                  % raw recording file path (empty if multiple files are sorted together)
    end

    %% LIFECYCLE
    methods
        function obj = Config(filename)
            %CONFIG Construct an instance of this class
            if nargin == 0
                return;
            end

            obj.isLoaded = false;
            obj.isError = false;

            if ~isfile(filename)
                emsg = sprintf('Cannot load config: file %s not found', filename);
                errordlg(emsg, 'Missing config file');
                obj.isError = true;
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

            s = jrclust.utils.mToStruct(obj.configFile);
            if isfield(s, 'template_file') && ~isempty(s.template_file)
                try
                    t = jrclust.utils.mToStruct(jrclust.utils.absPath(s.template_file));
                    fns = fieldnames(s);
                    for i = 1:numel(fns) % merge s (specific) into t (general)
                        t.(fns{i}) = s.(fns{i});
                    end
                    s = t;
                catch
                    wmsg = sprintf('Could not find template file %s', s.template_file);
                    warndlg(wmsg, 'Missing template file');
                end
            end

            fns = fieldnames(s);
            unusedProps = cell(numel(fns), 1);
            errorProps = cell(numel(fns), 2);

            % loop through all user-defined parameters
            for i = 1:numel(fns)
                % ignore configFile
                if strcmp(fns{i}, 'configFile') || strcmp(fns{i}, 'vcFile_prm')
                    continue;
                elseif strcmp(fns{i}, 'vcFile') % old single recording
                    if ~isempty(s.vcFile)
                        obj.setRawRecordings(s.vcFile);
                    end

                    continue;
                elseif strcmp(fns{i}, 'csFile_merge')
                    if ~isempty(s.csFile_merge)
                        obj.setRawRecordings(s.csFile_merge);
                    end

                    continue;
                end

                % empty values in the param file take on their defaults
                if ~isempty(s.(fns{i}))
                    if ~isprop(obj, fns{i})
                        unusedProps{i} = fns{i};
                        continue;
                    end

                    try
                        if ismember(fns{i}, {'probe_file', 'probeFile'})
                            s.(fns{i}) = jrclust.utils.absPath(s.(fns{i}), fileparts(obj.configFile));
                        end

                        obj.(fns{i}) = s.(fns{i});
                    catch ME % error or not a property we support
                        errorProps{i, 1} = fns{i};
                        errorProps{i, 2} = ME.message;
                    end
                end
            end

            % warn user of unrecognized props
            uu = cellfun(@(u) ~isempty(u), unusedProps);
            if any(uu)
                warnmsg = sprintf('The following properties were not recognized (possibly deprecated) and will be ignored:\n    %s', ...
                                  strjoin(unusedProps(uu), '\n    '));
                warndlg(warnmsg, 'Unrecognized properties');
            end

            ee = cellfun(@(e) ~isempty(e), errorProps(:, 1));
            if any(ee)
                emsgs = arrayfun(@(i) strjoin(errorProps(i, :), ': '), find(ee), 'UniformOutput', false);
                warnmsg = sprintf('These properties were not set due to errors:\n* %s', strjoin(emsgs, '\n\n * '));
                warndlg(warnmsg, 'Unset properties');
            end

            % check that exactly one of singleRaw and multiRaw is nonempty
            if isempty(obj.rawRecordings)
                errordlg('Specify rawRecordings (vcFile or csFile_merge)', 'Bad recording files');
                obj.isError = true;
                return;
            end

            % load probe from file
            if isempty(obj.probeFile)
                errordlg('Specify a probe file', 'Missing probe file');
                obj.isError = true;
                return;
            end
            try
                obj.loadProbe();
            catch ME
                errordlg(ME.message, 'Error loading probe file');
                obj.isError = true;
                return;
            end
        end

        function loadProbe(obj)
            %LOADPROBE Load probe construct from .mat or .m
            try % first try to load probe from mat file
                pstr = load(obj.probeFile, '-mat');
            catch
                pstr = jrclust.utils.mToStruct(obj.probeFile);
            end

            assert(isfield(pstr, 'channels'), 'missing `channels` field');
            obj.siteMap = pstr.channels;            

            assert(isfield(pstr, 'geometry'), 'missing `geometry` field');
            obj.siteLoc = pstr.geometry;

            assert(isfield(pstr, 'pad'), 'missing `pad` field');
            obj.probePad = pstr.pad;

            if ~isfield(pstr, 'shank')
                pstr.shank = [];
            end
            if ~isfield(pstr, 'cviShank')
                pstr.cviShank = [];
            end

            if isempty(pstr.shank) && ~isempty(pstr.cviShank)
                pstr.shank = pstr.cviShank;
            end
            if isempty(pstr.shank)
                obj.shankMap = ones(size(pstr.channels));
            elseif iscell(pstr.shank)
                obj.shankMap = cell2mat(arrayfun(@(i) i*ones(size(pstr.shank{i})), 1:numel(shank), 'UniformOutput', false));
            else
                obj.shankMap = pstr.shank;
            end

            if ~isempty(obj.nChans)
                obj.auxSites = setdiff(1:obj.nChans, 1:max(obj.siteMap));
            else
                obj.auxSites = [];
            end
        end

        function setRawRecordings(obj, rr)
            if ischar(rr)
                if ~isempty(dir(rr)) % first try current directory
                    dmr = {dir(rr)};
                    dmrDir = {dmr{1}.folder};
                    dmrName = {dmr{1}.name};
                    obj.rawRecordings = arrayfun(@(i) fullfile(dmrDir{i}, dmrName{i}), 1:numel(dmr{1}), 'UniformOutput', false);
                elseif ~isempty(dir(fullfile(fileparts(obj.configFile), rr))) % try relative to configFile directory
                    dmr = dir(fullfile(fileparts(obj.configFile), rr));
                    dmrDir = {dmr.folder};
                    dmrName = {dmr.name};
                    obj.rawRecordings = arrayfun(@(i) fullfile(dmrDir{i}, dmrName{i}), 1:numel(dmr), 'UniformOutput', false);
                else
                    error('could not find file(s) ''%''', rr);
                end
            elseif iscell(rr)
                % first try current directory
                filesExist = cellfun(@(m) isfile(jrclust.utils.absPath(m)), rr);
                if all(filesExist)
                    obj.rawRecordings = rr;
                    return;
                elseif any(filesExist) && ~all(filesExist) % some found but not all
                    errmsg = sprintf('The following recording files were specified but not found:\n    %s', ...
                                  strjoin(rr(filesExist), '\n    '));
                    errdlg(errmsg, 'Missing files');
                    error(errmsg);
                end

                % we have either returned or errored if we found any files
                % if we're still here, we have to look for them relative to
                % configFile's directory
                filesExist = cellfun(@(m) isfile(jrclust.utils.absPath(m, fileparts(obj.configFile))), rr);
                if all(filesExist)
                    obj.rawRecordings = rr;
                    return;
                elseif any(filesExist) && ~all(filesExist) % some found but not all
                    errmsg = sprintf('The following recording files were specified but not found:\n    %s', ...
                                  strjoin(rr(filesExist), '\n    '));
                    errdlg(errmsg, 'Missing files');
                    error(errmsg);
                end
            end

            % natsort filenames
            obj.rawRecordings = jrclust.utils.sortNat(obj.rawRecordings);
        end

        function validateParams(obj)
            %VALIDATEPARAMS Ensure parameters make sense
            % check probe is consistent
            if ~(jrclust.utils.ismatrixnum(obj.siteMap) && all(obj.siteMap > 0))
                errordlg('Malformed channel map', 'Bad probe configuration');
                obj.isError = true;
                return;
            end

            nc = numel(obj.siteMap);
            if ~(ismatrix(obj.siteLoc) && size(obj.siteLoc, 1) == nc)
                errordlg('Malformed probe geometry', 'Bad probe configuration');
                obj.isError = true;
                return;
            end

            if numel(obj.shankMap) ~= nc
                errordlg('Malformed shank indexing', 'Bad probe configuration');
                obj.isError = true;
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
                errordlg('nSitesExcl is too large or nSiteDir is too small', 'Bad configuration');
                obj.isError = true;
            end

            % try to infer a ground-truth file
            if isempty(obj.gtFile)
                try
                    obj.gtFile = strrep(obj.configFile, '.prm', '_gt.mat');
                catch % does not exist, leave empty
                end
            end

            % compute raw window limits if not given
            if isempty(obj.evtWindowRawms)
                % set in units of samples (will set ms units automatically)
                obj.evtWindowRawSamp = obj.evtWindowRawFactor * obj.evtWindowSamp;
            end
            
            obj.siteNeighbors = jrclust.utils.findSiteNeighbors(obj.siteLoc, 2*obj.nSiteDir + 1, obj.ignoreSites, obj.shankMap);

            % boost that gain
            obj.bitScaling = obj.bitScaling/obj.gainBoost;
        end
    end

    %% USER METHODS
    methods
        function val = getOr(obj, fn, dv)
            %GETOR GET set value obj.(fn) OR default value dv if unset or empty
            if nargin < 3
                dv = [];
            end

            if ~isprop(obj, fn) || isempty(obj.(fn))
                val = dv;
            else
                val = obj.(fn);
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

        % auxSites/viChan_aux
        function set.auxSites(obj, ac)
            assert(jrclust.utils.ismatrixnum(ac) && all(ac > 0), 'malformed auxSites');
            obj.auxSites = ac;
        end
        function ac = get.viChan_aux(obj)
            obj.logOldP('viChan_aux');
            ac = obj.auxSites;
        end
        function set.viChan_aux(obj, ac)
            obj.logOldP('viChan_aux');
            obj.auxSites = ac;
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
            bp = jrclust.utils.typeBytes(obj.dtype);
        end

        % carMode/vcCommonRef
        function set.carMode(obj, cm)
            legalTypes = {'mean', 'median', 'whiten', 'none'};
            failMsg = sprintf('legal carModes are %s', strjoin(legalTypes, ', '));
            assert(sum(strcmp(cm, legalTypes)) == 1, failMsg);
            obj.carMode = cm;
        end
        function cm = get.vcCommonRef(obj)
            obj.logOldP('vcCommonRef');
            cm = obj.carMode;
        end
        function set.vcCommonRef(obj, cm)
            obj.logOldP('vcCommonRef');
            obj.carMode = cm;
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
            legalTypes = {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'fftdiff', 'none'};
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

        % dtype/vcDataType
        function set.dtype(obj, dt)
            legalTypes = {'int16', 'uint16', 'int32', 'uint32', 'single', 'double'};
            failMsg = sprintf('legal dtypes are: %s', strjoin(legalTypes, ', '));
            assert(ismember(dt, legalTypes), failMsg);
            obj.dtype = dt;
        end
        function dt = get.vcDataType(obj)
            obj.logOldP('vcDataType');
            dt = obj.dtype;
        end
        function set.vcDataType(obj, dt)
            obj.logOldP('vcDataType');
            obj.dtype = dt;
        end

        % evtDetectRad/maxDist_site_um
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

        % evtManualThresh/spkThresh
        function mt = get.evtManualThresh(obj)
            mt = obj.evtManualThreshuV / obj.bitScaling;
        end
        function mt = get.spkThresh(obj)
            obj.logOldP('spkThresh');
            mt = obj.evtManualThresh;
        end

        % evtManualThreshuV/spkThresh_uV
        function set.evtManualThreshuV(obj, mt)
            assert(isempty(mt) || (jrclust.utils.isscalarnum(mt) && mt ~= 0)); % TODO: positive or negative?
            obj.evtManualThreshuV = mt;
        end
        function mt = get.spkThresh_uV(obj)
            obj.logOldP('spkThresh_uV');
            mt = obj.evtManualThreshuV;
        end
        function set.spkThresh_uV(obj, mt)
            obj.logOldP('spkThresh_uV');
            obj.evtManualThreshuV = mt;
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

        % evtWindowms/spkLim_ms
        function set.evtWindowms(obj, ew)
            assert(ismatrix(ew) && all(size(ew) == [1 2]) && ew(1) < 0 && ew(2) > 0, 'degenerate evtWindowms');
            obj.evtWindowms = ew;
        end
        function ew = get.spkLim_ms(obj)
            obj.logOldP('spkLim_ms');
            ew = obj.evtWindowms;
        end
        function set.spkLim_ms(obj, ew)
            obj.logOldP('spkLim_ms');
            obj.evtWindowms = ew;
        end

        % evtWindowRawFactor/spkLim_raw_factor
        function set.evtWindowRawFactor(obj, ef)
            assert(jrclust.utils.isscalarnum(ef) && ef > 0, 'evtWindowRawFactor must be a positive scalar');
            obj.evtWindowRawFactor = ef;
        end
        function ef = get.spkLim_raw_factor(obj)
            obj.logOldP('spkLim_raw_factor');
            ef = obj.evtWindowRawFactor;
        end
        function set.spkLim_raw_factor(obj, ef)
            obj.logOldP('spkLim_raw_factor');
            obj.evtWindowRawFactor = ef;
        end

        % evtWindowRawms/spkLim_raw_ms
        function set.evtWindowRawms(obj, ew)
            assert(ismatrix(ew) && all(size(ew) == [1 2]) && ew(1) < 0 && ew(2) > 0, 'degenerate evtWindowRawms');
            obj.evtWindowRawms = ew;
        end
        function ew = get.spkLim_raw_ms(obj)
            obj.logOldP('spkLim_raw_ms');
            ew = obj.evtWindowRawms;
        end
        function set.spkLim_raw_ms(obj, ew)
            obj.logOldP('spkLim_raw_ms');
            obj.evtWindowRawms = ew;
        end

        % evtWindowRawSamp/spkLim_raw
        function ew = get.evtWindowRawSamp(obj)
            ew = round(obj.evtWindowRawms * obj.sampleRate / 1000);
        end
        function set.evtWindowRawSamp(obj, ew)
            obj.evtWindowRawms = ew * 1000 / obj.sampleRate;
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
            ew = round(obj.evtWindowms * obj.sampleRate / 1000);
        end
        function set.evtWindowSamp(obj, ew)
            obj.evtWindowms = ew * 1000 / obj.sampleRate;
        end
        function ew = get.spkLim(obj)
            obj.logOldP('spkLim');
            ew = obj.evtWindowSamp;
        end
        function set.spkLim(obj, ew)
            obj.logOldP('spkLim');
            obj.evtWindowSamp = ew;
        end

        % fftThreshMAD/fft_thresh
        function set.fftThreshMAD(obj, ft)
            assert(jrclust.utils.isscalarnum(ft) && ft >= 0, 'fftThreshMAD must be a nonnegative scalar');
            obj.fftThreshMAD = ft;
        end
        function ft = get.fft_thresh(obj)
            obj.logOldP('fft_thresh');
            ft = obj.fftThreshMAD;
        end
        function set.fft_thresh(obj, ft)
            obj.logOldP('fft_thresh');
            obj.fftThreshMAD = ft;
        end

        % filtOrder
        function set.filtOrder(obj, fo)
            assert(jrclust.utils.isscalarnum(fo) && fo > 0, 'bad filtOrder');
            obj.filtOrder = fo;
        end

        % filterType/vcFilter
        function set.filterType(obj, ft)
            legalTypes = {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'fftdiff', 'none'};
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

        % fImportKsort/fImportKilosort
        function set.fImportKsort(obj, fi)
            obj.fImportKsort = true && fi;
        end
        function fi = get.fImportKilosort(obj)
            obj.logOldP('fImportKilosort');
            fi = obj.fImportKsort;
        end
        function set.fImportKilosort(obj, fi)
            obj.logOldP('fImportKilosort');
            obj.fImportKilosort = fi;
        end

        % freqLim
        function set.freqLim(obj, fl)
            assert(jrclust.utils.ismatrixnum(fl) && all(size(fl) == [1 2]) && all(fl >= 0), 'bad freqLim');
            obj.freqLim = fl;
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
            assert(jrclust.utils.isscalarnum(mb) && mb > 0, 'maxBytesLoad must be a positive scalar');
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

        % multiRaw/csFile_merge
        function mr = get.csFile_merge(obj)
            obj.logOldP('csFile_merge');
            mr = obj.multiRaw;
        end

        % nFet_use
        function set.nFet_use(obj, nf)
            assert(jrclust.utils.isscalarnum(nf) && ~isempty(intersect(nf, [1 2 3])), 'nFet_use must be 1, 2, or 3');
            obj.nFet_use = nf;
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
        function set.outputDir(obj, od)
            od_ = jrclust.utils.absPath(od);
            if isfile(od_)
                error('''%s'' is a file', od);
            end
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

        % refracIntms/spkRefrac_ms
        function set.refracIntms(obj, ri)
            assert(jrclust.utils.isscalarnum(ri) && ri > 0, 'refractory interval must be a positive scalar');
            obj.refracIntms = ri;
        end
        function ri = get.spkRefrac_ms(obj)
            obj.logOldP('spkRefrac_ms');
            ri = obj.refracIntms;
        end
        function set.spkRefrac_ms(obj, ri)
            obj.logOldP('spkRefrac_ms');
            obj.refracIntms = ri;
        end

        % refracIntSamp/spkRefrac
        function ri = get.refracIntSamp(obj)
            ri = round(obj.refracIntms * obj.sampleRate / 1000);
        end
        function set.refracIntSamp(obj, ri)
            obj.refracIntms = ri * 1000 / obj.sampleRate;
        end
        function ri = get.spkRefrac(obj)
            obj.logOldP('spkRefrac');
            ri = obj.refracIntSamp;
        end
        function set.spkRefrac(obj, ri)
            obj.logOldP('spkRefrac');
            obj.refracIntSamp = ri;
        end

        % repeatLower/fRepeat_clu
        function set.repeatLower(obj, rl)
            rl = rl && true;
            obj.repeatLower = rl;
        end
        function rl = get.fRepeat_clu(obj)
            obj.logOldP('fRepeat_clu');
            rl = obj.repeatLower;
        end
        function set.fRepeat_clu(obj, rl)
            obj.logOldP('fRepeat_clu');
            obj.repeatLower = rl;
        end

        % rlDetrendMode/vcDetrend_postclu
        function set.rlDetrendMode(obj, dm)
            legalTypes = {'global', 'local', 'logz', 'hidehiko', 'none'};
            failMsg = sprintf('legal rlDetrendModes are %s', strjoin(legalTypes, ', '));
            assert(sum(strcmp(dm, legalTypes)) == 1, failMsg);
            obj.rlDetrendMode = dm;
        end
        function dm = get.vcDetrend_postclu(obj)
            obj.logOldP('vcDetrend_postclu');
            dm = obj.rlDetrendMode;
        end
        function set.vcDetrend_postclu(obj, dm)
            obj.logOldP('vcDetrend_postclu');
            obj.rlDetrendMode = dm;
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
            obj.showRaw = true && sr;
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
                tf_ = jrclust.utils.absPath(tf);
                assert(isfile(tf_), 'could not find threshold file ''%s''', tf);
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

        % useElliptic/fEllip
        function set.useElliptic(obj, ue)
            assert(isscalar(ue));
            ue = logical(ue);
            obj.useElliptic = ue;
        end
        function ue = get.fEllip(obj)
            obj.logOldP('fEllip');
            ue = obj.useElliptic;
        end
        function set.fEllip(obj, ue)
            obj.logOldP('fEllip');
            obj.useElliptic = ue;
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
        function set.verbose(obj, vb)
            obj.verbose = true && vb;
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
