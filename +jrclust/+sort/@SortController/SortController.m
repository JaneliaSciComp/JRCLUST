classdef SortController < handle
    %SORTCONTROLLER Handle for sorting spikes into clusters using Rodriguez-Laio
    properties (Access=private, Transient)
        rhoCK;          % CUDA kernel for rho computation
        deltaCK;        % CUDA kernel for delta computation

        hCfg;           % Config object
    end

    properties(SetAccess=private, Transient)
        errMsg;         % error message, if any
        isError;        % flag indicating an error occurred in sorting
    end

    %% LIFECYCLE
    methods
        function obj = SortController(hCfg)
            %SORTCONTROLLER Constructor
            obj.hCfg = hCfg;
            obj.isError = 0;
        end
    end
end % class
