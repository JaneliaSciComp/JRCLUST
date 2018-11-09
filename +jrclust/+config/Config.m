classdef Config < handle & dynamicprops
    %CONFIG JRCLUST session configuration
    % replacement for P struct

    properties (SetObservable, SetAccess=private, Hidden)
        isLoaded;
        isError;
        oldPcount;
    end

    % old-style params, will be deprecated after a grace period,
    properties (SetObservable, Dependent, Hidden, Transient)
        % cvrDepth_drift;           % doesn't appear to be used (was {})
        % maxSite_detect;           % doesn't appear to be used (was 1.5)
        % maxSite_dip;              % doesn't appear to be used (was [])
        % maxSite_fet;              % doesn't appear to be used (was [])
        % maxSite_merge;            % doesn't appear to be used (was [])
        % maxSite_pix;              % doesn't appear to be used (was 4.5)
        % maxSite_show;             % appears to be synonymous with nSiteDir/maxSite
        % maxSite_sort;             % doesn't appear to be used (was [])
        % rejectSpk_mean_thresh;    % appears to be synonymous with blank_thresh
        blank_thresh;               % => blankThresh
        csFile_merge;               % => multiRaw
        fGpu;                       % => useGPU
        gain_boost;                 % => gainBoost;
        maxDist_site_um             % => evtMergeRad
        maxDist_site_spk_um;        % => evtDetectRad
        maxSite;                    % => nSiteDir
        miSites;                    % => siteNeighbors
        mrSiteXY;                   % => siteLoc
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
        tlim;                       % => dispTimeLimits
        tlim_load;                  % => loadTimeLimits
        uV_per_bit;                 % => bitScaling
        vcDataType;                 % => dtype
        vcFile;                     % => singleRaw
        vcFile_gt;                  % => gtFile
        vcFile_prm;                 % => configFile
        vcFilter;                   % => filterType
        vcFilter_show;              % => dispFilter
        viChan_aux;                 % => auxSites
        viShank_site;               % => shankMap
        vrSiteHW;                   % => probePad
        viSite2Chan;                % => siteMap
        viSiteZero;                 % => ignoreSites
    end

    % computed from other params but not in need of storage
    properties (SetObservable, Dependent)
        bytesPerSample;             % byte count for each discrete sample
        evtManualThresh;            % evtManualThreshuV / bitScaling
        evtWindowRawSamp;           % interval around event to extract raw spike waveforms, in samples
        evtWindowSamp;              % interval around event to extract filtered spike waveforms, in samples
        refracIntSamp;              % spike refractory interval, in samples
    end

    % new-style params
    properties (SetObservable)
        % computation params
        useGPU = true;              % use GPU in computation if true

        % recording params
        auxSites;                   %
        bitScaling = 0.30518;       % bit scaling factor (uV/bit)
        configFile;                 % parameter file
        dtype = 'int16';            % raw data binary format
        gainBoost = 1;              % multiply the raw recording by this gain to boost uV/bit
        gtFile= '';                 % ground truth file (default: SESSION_NAME_gt.mat) (TODO: specify format)
        lfpSampleRate = 2500;       % sample rate of the LFP recording, in Hz
        multiRaw;                   % list of recording files to merge (empty if single file is used)
        probeFile;                  % probe file to use (.prb, .mat)
        probePad;                   %
        sampleRate = 30000;         % sample rate of the recording, in Hz
        shankMap;                   % index of shank to which a site belongs
        singleRaw;                  % raw recording file path (empty if multiple files are sorted together)
        siteLoc;                    % x-y locations of channels on the probe, in microns
        siteMap;                    % channel mapping; row i in the data corresponds to channel `siteMap(i)`
        siteNeighbors;              % indices of neighbors for each site

        % preprocessing params
        loadTimeLimits = [];        % time range of recording to load, in s (use whole range if empty)
        filterType = 'ndiff';       % filter to use {'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'fftdiff', 'none'}

        % spike detection params
        blankThresh = [];           % reject spikes exceeding the channel mean after filtering (MAD unit), ignored if [] or 0
        evtDetectRad = 75;          % radius for extracting waveforms, in microns (used if nSiteDir and nSitesExcl are empty)
        evtManualThreshuV = [];     % manual spike detection threshold, in microvolts
        evtMergeRad = 50;           % radius of spike event merging, in microns
        evtWindowms = [-0.25 0.75]; % interval around event to extract filtered spike waveforms, in ms
        evtWindowRawms;             % interval around event to extract raw spike waveforms, in ms
        evtWindowRawFactor = 2;     % ratio of raw samples to filtered samples to extract if evtWindowRawms is not set
        ignoreSites = [];           % sites to manually ignore in the sorting
        nDiff_filt = 2;             % Differentiation filter for vcFilter='sgdiff', ignored otherwise. Set to [] to disable. 2n+1 samples used for centered differentiation
        nSiteDir;                   % number of neighboring sites to group in each direction (TODO: deprecate this)
        nSitesExcl;                 % number of sites to exclude from the spike waveform group
        refracIntms = 0.25;         % spike refractory interval, in ms
        siteCorrThresh = 0;         % reject bad sites based on max correlation with neighboring sites, using raw waveforms; ignored if 0

        % feature extraction params

        % clustering params
        log10RhoCut = -2.5;         % the base-10 log of the rho cutoff value

        % display params
        dispFilter = '';
        dispTimeLimits = [0 0.2];   % time range to display (in s?)

        % to get to, eventually
        LineStyle = '';
        MAX_BYTES_LOAD = [];
        MAX_LOAD_SEC = [];
        MAX_LOG = 5;
        S_imec3 = [];
        autoMergeCriterion = 'xcorr';
        blank_period_ms = 5;
        corrLim = [0.9 1];
        cviShank = [];
        dc_factor = 1;
        dc_frac = [];
        dc_percent = 2;
        delta1_cut = 0.6;
        dinput_imec_trial = 1;
        duration_file = [];
        fAddCommonRef = false;
        fAverageTrial_psth = true;
        fCacheRam = true;
        fCheckSites = false;
        fDetectBipolar = false;
        fDiscard_count = true;
        fDrift_merge = true;
        fEllip = true;
        fGroup_shank = false;
        fInterp_fet = true;
        fInverse_file = false;
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
        fRepeat_clu = false;
        fRms_detect = false;
        fRun = true;
        fSaveEvt = true;
        fSavePlot_RD = true;
        fSaveRawSpk = false;
        fSaveSpk = true;
        fShowAllSites = false;
        fSingleColumn_track = true;
        fSmooth_track = true;
        fSpatialMask_clu = false;
        fSpike_show = true;
        fText = true;
        fTranspose_bin = true;
        fUseCache_track = false;
        fUseLfp_track = true;
        fVerbose = false;
        fWav_raw_show = false;
        fWhiten_traces = false;
        fft_thresh = 0;
        filtOrder = 3;
        filter_sec_rate = 2;
        filter_shape_rate = 'triangle';
        flim_vid = [];
        freqLim = [300 3000];
        freqLimNotch = [];
        freqLimNotch_lfp = [];
        freqLim_corr = [15 150];
        freqLim_excl_track = [58 62];
        freqLim_lfp = [];
        freqLim_track = [15 150];
        header_offset = 0;
        iChan_aux = [];
        iChan_vid = [];
        iClu_show = [];
        iGpu = 1;
        load_fraction_track = [];
        maxAmp = 250;
        maxAmp_lfp = 1000;
        maxCluPerSite = 20;
        maxDist_site_merge_um = 35;
        maxLfpSdZ = 4.5;
        maxSite_track = [2 3 4 5 6 7 8];
        maxWavCor = 0.98;
        max_shift_track = [];
        min_count = 30;
        min_sites_mask = 5;
        mrColor_proj = [0.75 0.75 0.75; 0 0 0; 1 0 0];
        nBytes_file = [];
        nC_max = 45;
        nChans = 120;
        nClu_show_aux = 10;
        nInterp_merge = 1;
        nLoads_gpu = 8;
        nLoads_max_preview = 30;
        nMinAmp_ms = 0;
        nPad_filt = 100;
        nPcPerChan = 1;
        nPc_dip = 3;
        nRepeat_merge = 10;
        nShow = 200;
        nShow_proj = 500;
        nSites_excl_ref = 6;
        nSkip_show = 1;
        nSkip_whiten = 10;
        nSmooth_ms_psth = 50;
        nSpk_show = 30;
        nT_drift = [];
        nThreads = 128;
        nTime_clu = 4;
        nTime_traces = 1;
        nneigh_min_detect = 0;
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
        sort_file_merge = 1;
        spkLim_factor_merge = 1;
        spkLim_ms_fet = [-0.25 0.75];
        spkThresh_max_uV = [];
        tBin_track = 9;
        tRefrac_trial = 0.001;
        tbin_drift = [];
        tbin_psth = 0.01;
        template_file = '';
        thresh_automerge_pca = [];
        thresh_corr_track = [];
        thresh_mad_clu = 7.5;
        thresh_merge_clu = 0;
        thresh_sd_ref = 5;
        thresh_split_clu = 0;
        thresh_trial = [];
        time_tick_show = [];
        tlim_clu = [];
        tlim_lfp = [0 5];
        tlim_psth = [-1 5];
        tlim_vid = [];
        um_per_pix = 20;
        vcCluDist = 'eucldist';
        vcCluWavMode = 'mean';
        vcCommonRef = 'mean';
        vcDate_file = '';
        vcDc_clu = 'distr';
        vcDetrend_postclu = 'global';
        vcFet = 'gpca';
        vcFet_show = 'vpp';
        vcFile_aux = '';
        vcFile_bonsai = '';
        vcFile_lfp = '';
        vcFile_thresh = '';
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
        vnFilter_user = [];
        vrScale_aux = 1;
        xtick_psth = 0.2;
        ybin_drift = 2;
    end

    % lifecycle
    methods
        function obj = Config(filename)
            %CONFIG Construct an instance of this class
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

        function loadParams(obj)
            obj.oldPcount = containers.Map();

            s = jrclust.utils.mToStruct(obj.configFile);
            if ~isempty(s.template_file)
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

            % loop through all user-defined parameters
            for i = 1:numel(fns)
                % ignore configFile
                if strcmp(fns{i}, 'configFile') || strcmp(fns{i}, 'vcFile_prm')
                    continue;
                end

                % empty values in the param file take on their defaults
                if ~isempty(s.(fns{i}))
                    try
                        obj.(fns{i}) = s.(fns{i});
                    catch ME % error or not a property we support
                        unusedProps{i} = fns{i};
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

            % check that exactly one of singleRaw and multiRaw is nonempty
            if ~xor(isempty(obj.singleRaw), isempty(obj.multiRaw))
                errordlg('Specify exactly one of singleRaw (vcFile) or multiRaw (csFile_merge)', 'Bad recording files');
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

        function validateParams(obj)
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

    % getters/setters
    methods
        function logOldP(obj, prm)
            %LOGOLDP Increment the old-style parameter counter
            %   collect stats on usage of old parameters
            if ~isKey(obj.oldPcount, prm)
                obj.oldPcount(prm) = 1;
            else
                obj.oldPcount(prm) = obj.oldPcount(prm) + 1;
            end
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
            legalTypes = {'', 'ndiff', 'sgdiff', 'bandpass', 'fir1', 'user', 'fftdiff', 'none'};
            assert(sum(strcmp(ft, legalTypes)) == 1, 'legal filterTypes are: %s', strjoin(legalTypes, ', '));
            if isempty(ft)
                ft = obj.filterType;
            end
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
            assert(sum(strcmp(dt, legalTypes)) == 1, 'legal dtypes are: %s', strjoin(legalTypes, ', '));
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
            assert(jrclust.utils.isscalarnum(ed) && ed >= 0, 'evtDetectRad must be a nonnegative scalar');
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
            assert(jrclust.utils.isscalarnum(mt) && mt ~= 0); % TODO: positive or negative?
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
            ew = obj.evtWindowms;
        end
        function set.spkLim(obj, ew)
            obj.logOldP('spkLim');
            obj.evtWindowSamp = ew;
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
            gf_ = jrclust.utils.absPath(gf, fileparts(obj.configFile));
            assert(isfile(gf_), 'could not find recording file ''%s''', gf);
            obj.gtFile = gf_;
        end
        function gf = get.vcFile_gt(obj)
            obj.logOldP('vcFile_gt');
            gf = obj.gtFile;
        end
        function set.vcFile_gt(obj, gf)
            obj.logOldP('vcFile_gt');
            obj.gtFile = gf;
        end

        % ignoreSites/viSiteZero
        function set.ignoreSites(obj, ig)
            assert(ismatrix(ig) && all(ig > 0), 'degenerate ignoreSites');
            % don't manually ignore sites that are automatically ignored
            obj.ignoreSites = ig(ig > numel(obj.siteMap));
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

        % multiRaw/csFile_merge
        function set.multiRaw(obj, mr)
            if ischar(mr)
                if ~isempty(dir(mr)) % first try current directory
                    dmr = {dir(mr)};
                    dmrDir = {dmr{1}.folder};
                    dmrName = {dmr{1}.name};
                    obj.multiRaw = arrayfun(@(i) fullfile(dmrDir{i}, dmrName{i}), 1:numel(dmr{1}), 'UniformOutput', false);
                elseif ~isempty(dir(fullfile(obj.configFile, mr))) % try relative to configFile directory
                    dmr = dir(fullfile(obj.configFile, mr));
                    dmrDir = {dmr{1}.folder};
                    dmrName = {dmr{1}.name};
                    obj.multiRaw = arrayfun(@(i) fullfile(dmrDir{i}, dmrName{i}), 1:numel(dmr{1}), 'UniformOutput', false);
                else
                    error('could not find file(s) ''%''', mr);
                end
            elseif iscell(mr)
                % first try current directory
                filesExist = cellfun(@(m) isfile(jrclust.utils.absPath(m)), mr);
                if all(filesExist)
                    obj.multiRaw = mr;
                    return;
                elseif any(filesExist) && ~all(filesExist) % some found but not all
                    errmsg = sprintf('The following recording files were specified but not found:\n    %s', ...
                                  strjoin(mr(filesExist), '\n    '));
                    errdlg(errmsg, 'Missing files');
                    error(errmsg);
                end

                % we have either returned or errored if we found any files
                % if we're still here, we have to look for them relative to
                % configFile's directory
                filesExist = cellfun(@(m) isfile(jrclust.utils.absPath(m, fileparts(obj.configFile))), mr);
                if all(filesExist)
                    obj.multiRaw = mr;
                    return;
                elseif any(filesExist) && ~all(filesExist) % some found but not all
                    errmsg = sprintf('The following recording files were specified but not found:\n    %s', ...
                                  strjoin(mr(filesExist), '\n    '));
                    errdlg(errmsg, 'Missing files');
                    error(errmsg);
                end
            end
        end
        function mr = get.csFile_merge(obj)
            obj.logOldP('csFile_merge');
            mr = obj.multiRaw;
        end
        function set.csFile_merge(obj, mr)
            obj.logOldP('csFile_merge');
            obj.multiRaw = mr;
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

        % probeFile/probe_file
        function set.probeFile(obj, pf)
            pf_ = jrclust.utils.absPath(pf, fullfile(jrclust.utils.basedir, 'probes'));
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

        % shankMap/viSite2Chan
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

        % singleRaw/vcFile
        function set.singleRaw(obj, rf)
            rf_ = jrclust.utils.absPath(rf, fileparts(obj.configFile));
            assert(isfile(rf_), 'could not find recording file ''%s''', rf);
            obj.singleRaw = rf_;
        end
        function rf = get.vcFile(obj)
            obj.logOldP('vcFile');
            rf = obj.singleRaw;
        end
        function set.vcFile(obj, rf)
            obj.logOldP('vcFile');
            obj.singleRaw = rf;
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
            assert(jrclust.utils.ismatrixnum(sn) && all(size(sn) == [2*obj.nSiteDir + 1, numel(obj.siteMap)]), 'bad siteNeighbors');
            obj.siteNeighbors = sn;
        end
        function sn = get.miSites(obj)
            obj.logOldP('viSite2Chan');
            sn = obj.siteNeighbors;
        end
        function set.miSites(obj, sn)
            obj.logOldP('viSite2Chan');
            obj.siteNeighbors = sn;
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
    end
end

