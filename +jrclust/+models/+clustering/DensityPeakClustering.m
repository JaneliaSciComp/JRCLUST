classdef DensityPeakClustering < handle
    %DENSITYPEAKCLUSTERING Legacy class
    properties (Hidden, SetAccess=protected, SetObservable)
        sRes;               % sorting results
        dRes;               % detection results
    end

    properties (SetAccess=private, SetObservable)
        clusterCenters;     % cluster centers
        clusterCentroids;   % centroids of clusters on the probe
        unitCount;          % number of spikes per cluster
        clusterNotes;       % notes on clusters
        clusterSites;       % site on which spikes in this cluster most often occur
        editPos;            % current position in edit history
        history;            % cell array, log of merge/split/delete operations
        spikeClusters;      % individual spike assignments
        spikesByCluster;    % cell array of spike indices per cluster
        meanWfLocal;        % mean filtered waveforms for each cluster
        meanWfGlobal;       % mean filtered waveforms for each cluster over all sites
        meanWfLocalRaw;     % mean raw waveforms for each cluster
        meanWfGlobalRaw;    % mean raw waveforms for each cluster over all sites
        meanWfRawLow;       % mean raw waveforms for each cluster over all sites at a low point on the probe (for drift correction)
        meanWfRawHigh;      % mean raw waveforms for each cluster over all sites at a high point on the probe (for drift correction)

        nSitesOverThresh;   % number of sites exceeding the detection threshold, per cluster
        siteRMS;            % site-wise 
        unitPeaks;          % minimum voltage of mean filtered waveforms at peak site, per cluster
        unitPeaksRaw;       % minimum voltage (uV) of mean raw waveforms at peak site, per cluster
        unitPeakSites;      % sites on which unitPeaks occur
        unitVpp;            % peak-to-peak voltage of filtered waveforms at peak site, per cluster
        unitVppRaw;         % peak-to-peak voltage of raw waveforms at peak site, per cluster
        unitISIRatio;       % inter-spike interval ratio #(ISI <= 2ms)/#(ISI <= 20ms), per cluster
        unitIsoDist;        % isolation distance
        unitLRatio;         % L-ratio
        unitSNR;            % signal-to-noise ratio at peak site (peak/RMS)
    end

    %% LIFECYCLE
    methods
        function obj = DensityPeakClustering(varargin)
            warning('This is a legacy class');
        end
    end
end
