classdef SortController < handle
    %SORTCONTROLLER Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess=private, SetObservable, Transient)
        hCfg;
        isError;
    end

    % LIFECYCLE
    methods
        function obj = SortController(hCfg)
            obj.hCfg = hCfg;

            obj.isError = false;
        end
    end
end

