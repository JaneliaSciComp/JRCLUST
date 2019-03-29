classdef (Abstract) Clustering < handle
    %CLUSTERING A clustering of spike data
    %% CONFIGURATION
    properties (Hidden, SetObservable)
        hCfg;               % Config object
    end

    %% CLASS INTROSPECTION
    properties (SetAccess=protected)
        unitFields;         % data fields related to clusters
    end

    %% DETECTION/SORTING RESULTS
    properties (Hidden, SetObservable)
        dRes;               % detection results
        sRes;               % sorting results
    end

    %% SORTING DATA (MUTABLE)
    properties (SetObservable)
        clusterCentroids;   % centroids of clusters on the probe
        clusterNotes;       % notes on clusters
        clusterSites;       % mode site per cluster
        editPos;            % current position in edit history
        history;            % cell array, log of merge/split/delete operations
        meanWfGlobal;       % mean filtered waveforms for each cluster over all sites
        meanWfGlobalRaw;    % mean raw waveforms for each cluster over all sites
        meanWfLocal;        % mean filtered waveforms for each cluster
        meanWfLocalRaw;     % mean raw waveforms for each cluster
        meanWfRawLow;       % mean raw waveforms for each cluster over all sites at a low point on the probe (for drift correction)
        meanWfRawHigh;      % mean raw waveforms for each cluster over all sites at a high point on the probe (for drift correction)
        waveformSim;        % cluster similarity scores
        spikeClusters;      % individual spike assignments
        spikesByCluster;    % cell array of spike indices per cluster
        unitCount;          % number of spikes per cluster
    end

    % computed from other values, but only on set
    properties (SetAccess=protected, Transient)
        nClusters;          % number of clusters
    end

    % computed from other values
    properties (Dependent, Transient)
        nEdits;             % number of edits made to initial clustering
    end

    %% QUALITY METRICS
    properties (SetObservable)
        unitPeaks;          % minimum voltage of mean filtered waveforms at peak site, per cluster
        unitPeaksRaw;       % minimum voltage (uV) of mean raw waveforms at peak site, per cluster
        unitPeakSites;      % sites on which unitPeaks occur
        unitVpp;            % peak-to-peak voltage of filtered waveforms at peak site, per cluster
        unitVppRaw;         % peak-to-peak voltage of raw waveforms at peak site, per cluster
        unitISIRatio;       % inter-spike interval ratio #(ISI <= 2ms)/#(ISI <= 20ms), per cluster
        unitIsoDist;        % isolation distance
        unitLRatio;         % L-ratio
    end

    %% SORTING RESULTS (IMMUTABLE)
    properties (Dependent, Transient)
        initialClustering;  % initial assignment of spikes to cluster
    end

    %% DETECTION RESULTS (IMMUTABLE)
    properties (Dependent, Transient)
        spikeAmps;          % amplitudes of detected spikes
        spikePositions;     % positions on the probe at which spikes are detected
        spikeSites;         % sites on which spikes occur
        spikeTimes;         % times at which spikes occurred

        spikesBySite;       % aggregate of spike indices by site

        spikesRaw;          % raw spike traces
        spikesFilt;         % filtered spike traces
        spikeFeatures;      % features which were clustered
    end

    %% CACHED VALUES
    properties (Transient)
        spikesFiltVolt;     % spikesFilt in units of microvolts
        spikesRawVolt;      % spikesRaw in units of microvolts
    end

    %% LIFECYCLE
    methods
        function obj = Clustering()
            %CLUSTERING Construct an instance of this class
            fid = fopen(fullfile(jrclust.utils.basedir(), 'json', 'Clustering.json'), 'r');
            obj.unitFields = jsondecode(fread(fid, inf, '*char')');
            fclose(fid);

            obj.history = cell(0, 4);
        end
    end

    %% ABSTRACT METHODS
    methods (Abstract)
        autoMerge(obj);
        editSeek(obj, seekTo);
        success = exportQualityScores(obj, zeroIndex, fGui);
        rmOutlierSpikes(obj);
    end

    methods (Abstract, Access=protected, Hidden)
        nMerged = mergeBySim(obj);
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        postOp(obj, updateMe);
        removeEmptyClusters(obj);
    end

    %% GETTERS/SETTERS
    methods
        % hCfg
        function set.hCfg(obj, hc)
            failMsg = 'hCfg must be an object of type jrclust.Config';
            assert(isa(hc, 'jrclust.Config'), failMsg);
            obj.hCfg = hc;
        end

        % initialClustering
        function ic = get.initialClustering(obj)
            if isfield(obj.sRes, 'spikeClusters')
                ic = obj.sRes.spikeClusters;
            else
                ic = [];
            end
        end
        function set.initialClustering(obj, val)
            obj.sRes.spikeClusters = val;
        end

        % nEdits
        function ne = get.nEdits(obj)
            ne = size(obj.history, 1) - 1;
        end

        % spikeAmps
        function sa = get.spikeAmps(obj)
            if isfield(obj.dRes, 'spikeAmps')
                sa = obj.dRes.spikeAmps;
            else
                sa = [];
            end
        end
        function set.spikeAmps(obj, val)
            obj.dRes.spikeAmps = val;
        end

        % spikeClusters
        function set.spikeClusters(obj, sc)
            obj.spikeClusters = sc;
            obj.nClusters = numel(unique(sc(sc > 0))); %#ok<MCSUP>
        end

        % spikeFeatures
        function set.spikeFeatures(obj, sf)
            obj.dRes.spikeFeatures = sf;
        end
        function sf = get.spikeFeatures(obj)
            if isfield(obj.dRes, 'spikeFeatures')
                sf = obj.dRes.spikeFeatures;
            else
                sf = [];
            end
        end

        % spikePositions
        function sf = get.spikePositions(obj)
            if isfield(obj.dRes, 'spikePositions')
                sf = obj.dRes.spikePositions;
            else
                sf = [];
            end
        end
        function set.spikePositions(obj, val)
            obj.dRes.spikePositions = val;
        end

        % spikesBySite
        function ss = get.spikesBySite(obj)
            if isfield(obj.dRes, 'spikesBySite')
                ss = obj.dRes.spikesBySite;
            else
                ss = [];
            end
        end
        function set.spikesBySite(obj, val)
            obj.dRes.spikesBySite = val;
        end

        % spikesFilt
        function set.spikesFilt(obj, sf)
            obj.dRes.spikesFilt = sf;
        end
        function sf = get.spikesFilt(obj)
            if isfield(obj.dRes, 'spikesFilt')
                sf = obj.dRes.spikesFilt;
            else
                sf = [];
            end
        end

        % spikeSites
        function ss = get.spikeSites(obj)
            if isfield(obj.dRes, 'spikeSites')
                ss = obj.dRes.spikeSites;
            else
                ss = [];
            end
        end
        function set.spikeSites(obj, val)
            obj.dRes.spikeSites = val;
        end

        % spikesRaw
        function set.spikesRaw(obj, sr)
            obj.dRes.spikesRaw = sr;
        end
        function sr = get.spikesRaw(obj)
            if isfield(obj.dRes, 'spikesRaw')
                sr = obj.dRes.spikesRaw;
            else
                sr = [];
            end
        end

        % spikeTimes
        function val = get.spikeTimes(obj)
            if isfield(obj.dRes, 'spikeTimes')
                val = obj.dRes.spikeTimes;
            else
                val = [];
            end
        end
        function set.spikeTimes(obj, val)
            obj.dRes.spikeTimes = val;
        end
    end
end
