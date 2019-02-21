classdef DetectController < handle
    %DETECTCONTROLLER
    %% CONFIGURATION
    properties (SetObservable, Transient)
        hCfg;           % Config object
    end

    properties (SetAccess=private)
        errMsg;
        isError;
    end

    properties (Access=private)
        importTimes; % spike times already known
        importSites; % spike center sites already known

        siteThresh;
        spikeTimes;
        spikeAmps;
        centerSites;

        spikesRaw;
        spikesFilt;
        spikeFeatures;
    end

    %% LIFECYCLE
    methods
        function obj = DetectController(hCfg, importTimes, importSites)
            obj.hCfg = hCfg;
            if nargin < 2
                importTimes = [];
            end
            if nargin < 3
                importSites = [];
            end

            obj.importTimes = importTimes(:);
            obj.importSites = importSites(:);

            obj.isError = 0;
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        [spikeWindowsOut, spikeWindows2Out] = cancelOverlap(obj, spikeWindows, spikeWindows2, spikeTimes, spikeSites, spikeSites2, siteThresh);
        [spikeWindows, spikeTimes] = CARRealign(obj, spikeWindows, samplesIn, spikeTimes, neighbors);
        [windows, timeRanges] = extractWindows(obj, samplesIn, spTimes, spSites, fRaw);
        [spikeSites2, spikeSites3] = findSecondaryPeaks(obj, spikeWindows, spikeSites);
        spikeWindows = samplesToWindows2(obj, samplesIn, spikeSites, spikeTimes);
    end
end
