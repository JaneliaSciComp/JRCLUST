classdef CurateController < handle
    %CURATECONTROLLER Interface for manually curating sorted clusters

    properties (Access=private, Hidden, SetObservable)
        cRes;           % curate results struct, returned at endSession
        hCfg;           % Config object
        hClust;         % Clustering object
    end

    properties (Access=private, Hidden, Transient)
        hFigs;          % containers.Map of Figure objects
        selClusters;    % selected clusters, in order of selection
    end

    %% LIFECYCLE
    methods
        function obj = CurateController(hCfg)
            obj.hCfg = hCfg;
            obj.hFigs = containers.Map();
        end

        function delete(obj)
            if ~isempty(obj.hFigs)
                obj.closeFigures();
            end
        end
    end

    %% KEYPRESS/MOUSECLICK METHODS
    methods (Access=protected, Hidden)
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function closeFigures(obj)
            %CLOSEFIGURES Close all open figures
            hKeys = keys(obj.hFigs);
            for iKey = 1:numel(hKeys)
                hFig = obj.hFigs(hKeys{iKey});
                hFig.close();
            end

            obj.hFigs = containers.Map();
        end

        function plotAllFigures(obj)
            if isempty(obj.hFigs)
                obj.spawnFigures();
            end

            hKeys = keys(obj.hFigs);
            if any(strcmp(hKeys, 'hFigRD'))
                obj.hFigs('hFigRD') = doPlotFigRD(obj.hFigs('hFigRD'), obj.hClust, obj.hCfg);
            end
        end

        function spawnFigures(obj)
            %SPAWNFIGURES Create new figures
            obj.hFigs = doSpawnFigures(obj.hCfg);
        end
    end

    %% USER METHODS
    methods
        function beginSession(obj, hClust)
            %BEGINSESSION Start curating clusters
            obj.cRes = struct('hClust', hClust);
            obj.plotAllFigures();
        end

        function res = endSession(obj)
            %ENDSESSION Finish curating and return results
            res = obj.cRes;
            obj.cRes = [];
            obj.closeFigures();
        end
    end

    %% GETTERS/SETTERS
    methods
        function hc = get.hClust(obj)
            if ~isempty(obj.cRes) && isfield(obj.cRes, 'hClust')
                hc = obj.cRes.hClust;
            else
                hc = [];
            end
        end
    end
end

