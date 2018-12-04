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
            obj.selClusters = [];
        end

        function delete(obj)
            if ~isempty(obj.hFigs)
                obj.closeFigures();
            end
        end
    end

    %% KEYPRESS/MOUSECLICK METHODS
    methods (Hidden)
        function keyPressFigSim(obj, ~, hEvent)
            switch hEvent.Key
                case 'm' % merge
                    disp('ui_merge_');

                case 's' % split
                    disp('auto_split_(1)');

                case 'k' % kilosort template similarity view
                    % if get_set_(S0.P, 'fImportKsort', 0)
                    %     plot_FigWavCor_(S0, 'simscore');
                    % end
                    disp('kilosort');

                case {'d', 'backspace', 'delete'} % delete
                    disp('ui_delete_(S0)');

                case 'w' % waveform correlation view
                    % if get_set_(S0.P, 'fImportKsort', 0)
                    %     plot_FigWavCor_(S0, 'wavecor');
                    % end
                    disp('waveform');
            end
        end

        function mouseClickFigSim(obj, xyPoint, clickType)
            xyPoint = floor(xyPoint);
            if strcmp(clickType, 'normal') % left click
                obj.selClusters = xyPoint(1); % first selected cluster is x position

                if diff(xyPoint) ~= 0
                    obj.selClusters(2) = xyPoint(2);
                end

                disp(obj.selClusters);

                % S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0);
                % S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
            end
        end
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

            % plot rho-delta figure
            if any(strcmp(hKeys, 'hFigRD'))
                obj.hFigs('hFigRD') = doPlotFigRD(obj.hFigs('hFigRD'), obj.hClust, obj.hCfg);
            end

            % plot similarity score figure
            if any(strcmp(hKeys, 'hFigSim'))
                hFigSim = doPlotFigSim(obj.hFigs('hFigSim'), obj.hClust, obj.hCfg);
                hFigSim.hFunKey = @obj.keyPressFigSim;
                hFigSim.setMouseable(@obj.mouseClickFigSim);
                obj.hFigs('hFigSim') = hFigSim;
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

