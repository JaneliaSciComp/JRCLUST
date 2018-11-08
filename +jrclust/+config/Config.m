classdef Config < handle & dynamicprops
    %CONFIG JRCLUST session configuration
    % replacement for P struct
    
    properties (SetObservable, SetAccess=private, Hidden)
        isCompleted;
        isError;
        depCtr;
    end

    % old-style params, will be deprecated after a grace period,
    properties (SetObservable, Dependent, Hidden, Transient)
        rho_cut; % log10RhoCut
        vcFile; % singleRaw
        vcFile_prm; % configFile
    end

    % new-style params
    properties (SetObservable)
        % files
        configFile;         % parameter file
        singleRaw;          % raw recording file path (empty if multiple files are sorted together)
        
        % clustering params
        log10RhoCut = -2.5; % the base-10 log of the rho cutoff value

        % to get to, eventually
        LineStyle = '';
        MAX_BYTES_LOAD = [];
        MAX_LOAD_SEC = [];
        MAX_LOG = 5;
        S_imec3 = [];
        autoMergeCriterion = 'xcorr';
        blank_period_ms = 5;
        blank_thresh = [];
        corrLim = [0.9 1];
        csFile_merge = [];
        cviShank = [];
        cvrDepth_drift = [];
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
        fGpu = true;
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
        gain_boost = 1;
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
        maxDist_site_spk_um = 75;
        maxDist_site_um = 50;
        maxLfpSdZ = 4.5;
        maxSite = [];
        maxSite_detect = 1.5;
        maxSite_dip = [];
        maxSite_fet = [];
        maxSite_merge = [];
        maxSite_pix = 4.5;
        maxSite_show = [];
        maxSite_sort = [];
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
        nDiff_filt = 2;
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
        nSites_ref = [];
        nSkip_lfp = [];
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
        probe_file = '';
        qqFactor = 5;
        qqSample = 4;
        rateLim_psth = [];
        refrac_factor = 2;
        rms_filt_ms = 0;
        sRateHz = 25000;
        sRateHz_aux = [];
        sRateHz_lfp = 2500;
        sRateHz_rate = 1000;
        sec_per_load_preview = 1;
        slopeLim_ms = [0.05 0.35];
        sort_file_merge = 1;
        spkLim_factor_merge = 1;
        spkLim_ms = [-0.25 0.75];
        spkLim_ms_fet = [-0.25 0.75];
        spkLim_raw_factor = 2;
        spkLim_raw_ms = [];
        spkRefrac_ms = 0.25;
        spkThresh_max_uV = [];
        spkThresh_uV = [];
        tBin_track = 9;
        tRefrac_trial = 0.001;
        tbin_drift = [];
        tbin_psth = 0.01;
        template_file = '';
        thresh_automerge_pca = [];
        thresh_corr_bad_site = [];
        thresh_corr_track = [];
        thresh_mad_clu = 7.5;
        thresh_merge_clu = 0;
        thresh_sd_ref = 5;
        thresh_split_clu = 0;
        thresh_trial = [];
        time_tick_show = [];
        tlim = [0 0.2];
        tlim_clu = [];
        tlim_lfp = [0 5];
        tlim_load = [];
        tlim_psth = [-1 5];
        tlim_vid = [];
        uV_per_bit = 0.30518;
        um_per_pix = 20;
        vcCluDist = 'eucldist';
        vcCluWavMode = 'mean';
        vcCommonRef = 'mean';
        vcDataType = 'int16';
        vcDate_file = '';
        vcDc_clu = 'distr';
        vcDetrend_postclu = 'global';
        vcFet = 'gpca';
        vcFet_show = 'vpp';
        vcFile_aux = '';
        vcFile_bonsai = '';
        vcFile_gt = '';
        vcFile_lfp = '';
        vcFile_thresh = '';
        vcFile_trial = '';
        vcFile_vid = '';
        vcFilter = 'ndiff';
        vcFilter_detect = '';
        vcFilter_show = '';
        vcLabel_aux = '';
        vcMode_track = 'mt_cpsd2_mr';
        vcSpatialFilter = 'none';
        vcSpkRef = 'nmean';
        viChan_aux = [];
        viChan_bin = [];
        viChan_show = [];
        viDepth_excl_track = [];
        viDepth_track = [];
        viSiteZero = [];
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
            obj.isCompleted = false;
            obj.isError = false;

            if ~isfile(filename)
                emsg = sprintf('cannot load config: file %s not found', filename);
                errordlg(emsg, 'Missing config file');
                return;
            elseif isempty(which(filename))
                obj.configFile = filename;
            else
                obj.configFile = which(filename);
            end

            obj.loadParams();
        end

        function loadParams(obj)
            obj.depCtr = containers.Map();

            s = jrclust.utils.mToStruct(obj.configFile);
            fns = fieldnames(s);
            
            unusedProps = cell(numel(fns), 1);

            for i = 1:numel(fns)
                % ignore configFile
                if strcmpi(fns{i}, 'configFile') || strcmpi(fns{i}, 'vcFile_prm')
                    continue;
                end

                try
                    obj.(fns{i}) = s.(fns{i});
                catch ME % error or not a property we support
                    unusedProps{i} = fns{i};
                end
            end

            % warn user of unrecognized props
            uu = cellfun(@(u) ~isempty(u), unusedProps);
            if any(uu)
                warnmsg = sprintf('The following properties were not recognized and will be ignored:\n    %s', ...
                                  strjoin(unusedProps(uu), '\n    '));
                warndlg(warnmsg, 'Unrecognized properties');
            end
        end
    end

    % getters/setters
    methods
        function incDepCtr(obj, prm)
            % having spotted a deprecated parameter, note it
            if ~isKey(obj.depCtr, prm)
                obj.depCtr(prm) = 1;
            else
                obj.depCtr(prm) = obj.depCtr(prm) + 1;
            end
        end

        % configFile/vcFile_prm
        function set.configFile(obj, cf)
            cf_ = jrclust.utils.absPath(cf);
            assert(isfile(cf_), 'could not find file ''%s''', cf);
            obj.configFile = cf_;
        end
        function cf = get.vcFile_prm(obj)
            obj.incDepCtr('vcFile_prm');
            cf = obj.configFile;
        end
        function set.vcFile_prm(obj, cf)
            obj.incDepCtr('vcFile_prm');
            obj.configFile = cf;
        end

        % log10RhoCut/rho_cut
        function set.log10RhoCut(obj, rc)
            assert(isnumeric(rc) && isscalar(rc), 'log10RhoCut must be a numeric scalar');
            obj.log10RhoCut = rc;
        end
        function rc = get.rho_cut(obj)
            obj.incDepCtr('rho_cut');
            rc = obj.log10RhoCut;
        end
        function set.rho_cut(obj, rc)
            obj.incDepCtr('rho_cut');
            obj.log10RhoCut = rc;
        end

        % singleRaw/vcFile
        function set.singleRaw(obj, rf)
            rf_ = jrclust.utils.absPath(rf, fileparts(obj.configFile));
            assert(isfile(rf_), 'could not find file ''%s''', rf);
            obj.singleRaw = rf_;
        end
        function rf = get.vcFile(obj)
            obj.incDepCtr('vcFile');
            rf = obj.singleRaw;
        end
        function set.vcFile(obj, rf)
            obj.incDepCtr('vcFile');
            obj.singleRaw = rf;
        end
    end
end

