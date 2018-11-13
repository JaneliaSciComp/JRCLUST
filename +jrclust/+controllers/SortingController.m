classdef SortingController < handle
    %SORTINGCONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private, SetObservable, Transient)
        hCfg;
    end
    
    methods
        function obj = SortingController(hCfg)
            obj.hCfg = hCfg;
        end
    end
end

