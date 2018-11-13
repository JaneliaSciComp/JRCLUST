classdef DetectionController
    %DETECTIONCONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        hCfg;
        hRecs;
        nRecs;
    end
    
    % LIFECYCLE
    methods
        function obj = DetectionController(hCfg)
            obj.hCfg = hCfg;
            obj.nRecs = numel(obj.hCfg.rawRecordings);
            obj.hRecs = jrclust.models.Recording.empty;
        end
    end

    % USER METHODS
    methods
        function res = detect(obj, givenTimes, givenSites)
            res = struct();
            t0 = tic;

            if nargin < 2
                givenTimes = [];
            end
            if nargin < 3
                givenSites = [];
            end

            givenTimes = givenTimes(:);
            givenSites = givenSites(:);

            for i = 1:obj.nRecs
                fn = obj.hCfg.rawRecordings{i};
                disp(fn);
            end

            res.runtime = toc(t0);
        end
    end

    % UTILITY METHODS
    methods (Access=protected, Hidden)
        function [nLoad1, nSamples_load1, nSamples_last1] = planLoad(obj, nBytes_file, P)
            % plan load file size according to the available memory and file size (nBytes_file1)
            LOAD_FACTOR = 5; %GPU memory usage factor. 4x means 1/4 of GPU memory can be loaded

            nSamples1 = floor(nBytes_file / bytesPerSample_(P.vcDataType) / P.nChans);
            % nSamples_max = floor(mem_max_(P) / P.nChans / 4); % Bound by MAX_BYTES_LOAD
            if ~isfield(P, 'MAX_BYTES_LOAD'), P.MAX_BYTES_LOAD = []; end
            if isempty(P.MAX_BYTES_LOAD), P.MAX_BYTES_LOAD = floor(mem_max_(P) / LOAD_FACTOR); end
            if isempty(P.MAX_LOAD_SEC)
                nSamples_max = floor(P.MAX_BYTES_LOAD / P.nChans / bytesPerSample_(P.vcDataType));
            else
                nSamples_max = floor(P.sRateHz * P.MAX_LOAD_SEC);
            end

            if ~P.fTranspose_bin %load all in one, Catalin's format
                [nLoad1, nSamples_load1, nSamples_last1] = deal(1, nSamples1, nSamples1);
            else
                [nLoad1, nSamples_load1, nSamples_last1] = partition_load_(nSamples1, nSamples_max);
            end
        end %func

    end
end

