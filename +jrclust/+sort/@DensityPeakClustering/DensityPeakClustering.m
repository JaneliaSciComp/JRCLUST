classdef DensityPeakClustering < jrclust.interfaces.Clustering
    %DENSITYPEAKCLUSTERING A Rodriguez-Laio clustering of spike data
    %% OLD-STYLE PROPERTIES, publicly gettable (will be deprecated after a grace period)
    properties (SetObservable, Dependent, Hidden, Transient)
        cviSpk_clu;         % => spikesByCluster
        csNote_clu;         % => clusterNotes
        delta;              % => spikeDelta
        icl;                % => clusterCenters
        mrWavCor;           % => waveformSim
        nClu;               % => nClusters
        nneigh;             % => spikeNeigh
        ordrho;             % => ordRho
        rho;                % => spikeRho
        tmrWav_clu;         % => meanWfGlobal
        tmrWav_spk_clu;     % => meanWfGlobal
        tmrWav_raw_clu;     % => meanWfGlobalRaw
        trWav_spk_clu;      % => meanWfLocal
        trWav_raw_clu;      % => meanWfLocalRaw
        tmrWav_raw_hi_clu;  % => meanWfRawHigh
        tmrWav_raw_lo_clu;  % => meanWfRawLow
        viClu;              % => spikeClusters
        viClu_auto;         % => initialClustering
        viSite_clu;         % => clusterSites
        viSite_min_clu;     % => unitPeakSites
        vnSite_clu;         % => nSitesOverThresh
        vnSpk_clu;          % => unitCount
        vrIsiRatio_clu;     % => unitISIRatio
        vrIsoDist_clu;      % => unitIsoDist
        vrLRatio_clu;       % => unitLRatio
        vrPosX_clu;         % => clusterCentroids(:, 1)
        vrPosY_clu;         % => clusterCentroids(:, 2)
        vrSnr_clu;          % => unitSNR
        vrVmin_clu;         % => unitPeaks
        vrVmin_uv_clu;      % => unitPeaksRaw
        vrVpp_clu;          % => unitVpp
        vrVpp_uv_clu;       % => unitVppRaw
        vrVrms_site;        % => siteRMS
    end

    %% NEW-STYLE CLUSTERING PROPERTIES
    properties (SetObservable)
        clusterCenters;     % cluster centers
    end

    %% QUALITY METRICS
    properties (SetObservable)
        nSitesOverThresh;   % number of sites exceeding the detection threshold, per cluster
        siteRMS;            % site-wise threshold/qqFactor
        unitSNR;            % signal-to-noise ratio at peak site (peak/RMS)
    end

    %% SORTING RESULTS (IMMUTABLE)
    properties (Dependent, Transient)
        ordRho;             % spike-wise index ordered by density (DPCLUS)
        rhoCuts;            % site-wise distance cutoff values for computing density (DPCLUS)
        spikeDelta;         % spike-wise distance to nearest neighbor of higher density (DPCLUS)
        spikeNeigh;         % nearest neighbor of higher density (DPCLUS)
        spikeRho;           % spike-wise density (DPCLUS)
    end

    %% DETECTION RESULTS (IMMUTABLE)
    properties (Dependent, Transient)
        meanSiteThresh;     % mean sitewise detection threshold over all chunks
        siteThresh;         % sitewise detection threshold over all chunks
        spikesBySite2;      % aggregate of secondary spike indices by site
        spikesBySite3;      % aggregate of tertiary spike indices by site
        spikeSites2;        % secondary sites on which spikes occur
    end

    %% LIFECYCLE
    methods
        function obj = DensityPeakClustering(sRes, dRes, hCfg)
            fid = fopen(fullfile(jrclust.utils.basedir(), 'json', 'DensityPeakClustering.json'), 'r');
            dpFields = jsondecode(fread(fid, inf, '*char')');
            fclose(fid);
            obj.unitFields.vectorFields = [obj.unitFields.vectorFields; dpFields.vectorFields];

            obj.dRes = dRes;
            obj.hCfg = hCfg;
            isImport = obj.tryImport(sRes);

            if ~isImport
                obj.sRes = sRes;
                obj.spikeClusters = obj.initialClustering;
            end

            % these fields are mutable so we need to store copies in obj
            if isfield(sRes, 'clusterCenters')
                obj.clusterCenters = sRes.clusterCenters;
            else
                obj.clusterCenters = [];
            end
            if isfield(sRes, 'clusterCentroids')
                obj.clusterCentroids = sRes.clusterCentroids;
            else
                obj.clusterCentroids = [];
            end

            if ~isImport
                obj.clearNotes();
                obj.refresh(1, []);
                commitMsg = sprintf('%s;initial commit', datestr(now, 31));
                obj.commit(commitMsg);
            end
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        [sites1, sites2, sites3] = getSecondaryPeaks(obj);
        nMerged = mergeBySim(obj);
        postOp(obj, updateMe);
        success = tryImport(obj, sRes);
    end

    %% GETTERS/SETTERS
    methods
        % clusterCenters/icl
        function vals = get.icl(obj)
            vals = obj.clusterCenters;
        end

        % clusterNotes/csNote_clu
        function cn = get.csNote_clu(obj)
            cn = obj.clusterNotes;
        end

        % unitCount/vnSpk_clu
        function cc = get.vnSpk_clu(obj)
            cc = obj.unitCount;
        end

        % clusterSites/viSite_clu
        function cs = get.viSite_clu(obj)
            cs = obj.clusterSites;
        end

        % initialClustering/viClu_auto
        function ic = get.viClu_auto(obj)
            ic = obj.initialClustering;
        end

        % meanWfGlobal/tmrWav_spk_clu
        function mw = get.tmrWav_spk_clu(obj)
            mw = obj.meanWfGlobal;
        end

        % meanWfGlobalRaw/tmrWav_raw_clu
        function mw = get.tmrWav_raw_clu(obj)
            mw = obj.meanWfGlobalRaw;
        end

        % meanWfLocal/trWav_spk_clu
        function mw = get.trWav_spk_clu(obj)
            mw = obj.meanWfLocal;
        end

        % meanWfLocalRaw/trWav_raw_clu
        function mw = get.trWav_raw_clu(obj)
            mw = obj.meanWfLocalRaw;
        end

        % meanWfRawHigh/tmrWav_raw_hi_clu
        function mw = get.tmrWav_raw_hi_clu(obj)
            mw = obj.meanWfRawHigh;
        end

        % meanWfRawLow/tmrWav_raw_lo_clu
        function mw = get.tmrWav_raw_lo_clu(obj)
            mw = obj.meanWfRawLow;
        end

        % unitPeakSites/viSite_min_clu
        function ms = get.viSite_min_clu(obj)
            ms = obj.unitPeakSites;
        end

        % meanSiteThresh
        function st = get.meanSiteThresh(obj)
            if isfield(obj.dRes, 'meanSiteThresh')
                st = obj.dRes.meanSiteThresh;
            else
                st = [];
            end
        end
        function set.meanSiteThresh(obj, val)
            obj.dRes.meanSiteThresh = val;
        end

        % nClusters/nClu
        function nc = get.nClu(obj)
            nc = obj.nClusters;
        end

        % ordRho/ordrho
        function val = get.ordRho(obj)
            if isfield(obj.sRes, 'ordRho')
                val = obj.sRes.ordRho;
            else
                val = [];
            end
        end
        function val = get.ordrho(obj)
            val = obj.ordRho;
        end
        function set.ordRho(obj, val)
            obj.sRes.ordRho = val;
        end

        % rhoCuts
        function val = get.rhoCuts(obj)
            if isfield(obj.sRes, 'rhoCutSite')
                val = obj.sRes.rhoCutSite;
            else
                val = [];
            end
        end
        function set.rhoCuts(obj, val)
            obj.sRes.rhoCutSite = val;
        end

        % waveformSim/mrWavCor
        function ss = get.mrWavCor(obj)
            ss = obj.waveformSim;
        end

        % nSitesOverThresh/vnSite_clu
        function so = get.vnSite_clu(obj)
            so = obj.nSitesOverThresh;
        end

        % siteThresh
        function st = get.siteThresh(obj)
            if isfield(obj.dRes, 'meanSiteThresh')
                st = obj.dRes.meanSiteThresh;
            elseif isfield(obj.dRes, 'siteThresh') % backwards compatibility
                st = obj.dRes.siteThresh;
            else
                st = [];
            end
        end
        function set.siteThresh(obj, val)
            obj.dRes.siteThresh = val;
        end

        % spikeClusters/viClu
        function sc = get.viClu(obj)
            sc = obj.spikeClusters;
        end

        % spikeDelta/delta
        function sd = get.spikeDelta(obj)
            if isfield(obj.sRes, 'spikeDelta')
                sd = obj.sRes.spikeDelta;
            else
                sd = [];
            end
        end
        function sd = get.delta(obj)
            sd = obj.spikeDelta;
        end
        function set.spikeDelta(obj, val)
            obj.sRes.spikeDelta = val;
        end

        % spikeNeigh/nneigh
        function sn = get.spikeNeigh(obj)
            if isfield(obj.sRes, 'spikeNeigh')
                sn = obj.sRes.spikeNeigh;
            else
                sn = [];
            end
        end
        function sd = get.nneigh(obj)
            sd = obj.spikeNeigh;
        end
        function set.spikeNeigh(obj, val)
            obj.sRes.spikeNeigh = val;
        end

        % spikeRho/rho
        function sr = get.spikeRho(obj)
            if isfield(obj.sRes, 'spikeRho')
                sr = obj.sRes.spikeRho;
            else
                sr = [];
            end
        end
        function sr = get.rho(obj)
            sr = obj.spikeRho;
        end
        function set.spikeRho(obj, val)
            obj.sRes.spikeRho = val;
        end

        % spikesBySite2
        function ss = get.spikesBySite2(obj)
            if isfield(obj.dRes, 'spikesBySite2')
                ss = obj.dRes.spikesBySite2;
            else
                ss = [];
            end
        end
        function set.spikesBySite2(obj, val)
            obj.dRes.spikesBySite2 = val;
        end

        % spikesBySite3
        function ss = get.spikesBySite3(obj)
            if isfield(obj.dRes, 'spikesBySite3')
                ss = obj.dRes.spikesBySite3;
            else
                ss = [];
            end
        end
        function set.spikesBySite3(obj, val)
            obj.dRes.spikesBySite3 = val;
        end

        % spikeSites2
        function ss = get.spikeSites2(obj)
            if isfield(obj.dRes, 'spikeSites2')
                ss = obj.dRes.spikeSites2;
            else
                ss = [];
            end
        end
        function set.spikeSites2(obj, val)
            obj.dRes.spikeSites2 = val;
        end

        % tmrWav_clu
        function tm = get.tmrWav_clu(obj)
            tm = obj.meanWfGlobal;
        end

        % unitISIRato/VrIsiRatio_clu
        function ir = get.vrIsiRatio_clu(obj)
            ir = obj.unitISIRatio;
        end

        % unitIsoDist/vrIsoDist_clu
        function id = get.vrIsoDist_clu(obj)
            id = obj.unitIsoDist;
        end

        % unitLRatio/vrLRatio_clu
        function lr = get.vrLRatio_clu(obj)
            lr = obj.unitLRatio;
        end

        % unitPeaks/vrVmin_clu
        function pv = get.vrVmin_clu(obj)
            pv = obj.unitPeaks;
        end

        % unitPeaksRaw/vrVmin_uv_clu
        function pv = get.vrVmin_uv_clu(obj)
            pv = obj.unitPeaksRaw;
        end

        % unitSNR/vrSnr_clu
        function sn = get.vrSnr_clu(obj)
            sn = obj.unitSNR;
        end

        % unitVpp/vrVpp_clu
        function pv = get.vrVpp_clu(obj)
            pv = obj.unitVpp;
        end

        % unitVppRaw/vrVpp_uv_clu
        function pv = get.vrVpp_uv_clu(obj)
            pv = obj.unitVppRaw;
        end

        % vrPos{X,Y}_clu
        function xp = get.vrPosX_clu(obj)
            if ~isempty(obj.clusterCentroids)
                xp = obj.clusterCentroids(:, 1);
            else
                xp = [];
            end
        end
        function yp = get.vrPosY_clu(obj)
            if ~isempty(obj.clusterCentroids)
                yp = obj.clusterCentroids(:, 2);
            else
                yp = [];
            end
        end
    end
end
