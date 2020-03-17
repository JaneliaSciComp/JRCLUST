classdef CurateController < handle
    %CURATECONTROLLER Interface for manually curating sorted clusters
    %% CONFIGURATION
    properties (Dependent)
        hCfg;           % Config object
        hClust;         % Clustering object
        nShown;         % number of units actually shown
    end

    properties (SetAccess=private, Hidden, SetObservable)
        res;            % detect/sort results struct, passed by JRC controller
        cRes;           % curate results struct, returned at endSession
    end

    properties (AbortSet, SetAccess=private, Hidden, Transient, SetObservable)
        currentSite;    % current site in FigTime
        helpTexts;      % figure help texts
        hFigs;          % containers.Map of Figure objects
        hMenus;         % containers.Map of Menu handles
        isEnding;       % to prevent double prompting when endSession is called
        isWorking;      % to prevent keypress/mouseclick functions from colliding
        maxAmp;         % scaling factor for
        projSites;      % current sites in FigProj
        showSubset;     % subset of units to display
    end

    properties (SetAccess=private, Hidden, Transient, SetObservable)
        selected;       % selected clusters, in order of selection
    end

    properties (SetAccess=private, Hidden, Transient, SetObservable)
        spatial_idx;      % idx to sort the channels by their desired plotting order
        channel_idx;      % idx to sort back to channel number (useful internally)
    end

    %% LIFECYCLE
    methods
        function obj = CurateController(res)
            % transfer hClust from res to cRes
            obj.hClust = res.hClust;
            obj.res = rmfield(res, 'hClust');

            obj.hFigs = containers.Map();
            obj.hMenus = containers.Map();
            obj.isEnding = 0;
            obj.isWorking = 0;
            obj.currentSite = [];
            obj.maxAmp = [];
            obj.selected = [];
            obj.showSubset = 1:obj.hClust.nClusters;

            % load helptexts
            helpFile = fullfile(jrclust.utils.basedir(), 'json', 'helptexts.json');
            if exist(helpFile, 'file') == 2
                fid = fopen(helpFile, 'r');
                obj.helpTexts = jsondecode(fread(fid, inf, '*char')');
                fclose(fid);
            end
            obj.spatial_idx= 1:length(obj.hClust.spikesBySite);
        end

        function delete(obj)
            if ~isempty(obj.hFigs)
                obj.closeFigures();
            end
        end
    end

    %% KEYPRESS/MOUSECLICK METHODS
    methods (Hidden)
        keyPressFigSim(obj, hObject, hEvent);
        keyPressFigProj(obj, hObject, hEvent);
        keyPressFigTime(obj, hObject, hEvent);
        keyPressFigWav(obj, hObject, hEvent);
        mouseClickFigSim(obj, xyPos, clickType);
        mouseClickFigWav(obj, xyPos, clickType);
    end

    %% SPECIFIC FIGURE METHODS
    methods (Hidden)
        updateCursorFigSim(obj);
        updateCursorFigWav(obj);
        updateFigCorr(obj);
        updateFigHist(obj);
        updateFigISI(obj);
        updateFigMap(obj);
        updateFigPos(obj);
        updateFigProj(obj, doAutoscale);
        updateFigRD(obj);
        updateFigTime(obj, doAutoscale);
        updateFigPSTH(obj, doCreate);
        updateFigWav(obj);
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        addMenu(obj, hFig);
        closeFigures(obj);
        deleteAnnotated(obj);
        res = figApply(obj, hFun, varargin);
        killFigWav(obj, hObject, hEvent);
        plotAllFigures(obj);
        spawnFigures(obj);
        updateHistMenu(obj);
        updateNoteMenu(obj);
    end

    methods (Static)
        figPos = getCurrentFigPos();
        setFigPos(figPos);
    end


    %% GETTERS/SETTERS
    methods
        % hCfg
        function val = get.hCfg(obj)
            if ~isempty(obj.hClust)
                val = obj.hClust.hCfg;
            else
                val = [];
            end
        end

        % hClust
        function hc = get.hClust(obj)
            if ~isempty(obj.cRes) && isfield(obj.cRes, 'hClust')
                hc = obj.cRes.hClust;
            else
                hc = [];
            end
        end
        function set.hClust(obj, val)
            obj.cRes.hClust = val;
        end

        function set.spatial_idx(obj, val)
            if isempty(obj.spatial_idx)
                obj.spatial_idx = val;
                [~,obj.channel_idx] = sort(obj.spatial_idx);                
            else
                nChanged = sum(obj.spatial_idx(:)~=val(:));
                if ~nChanged
                    jrclust.utils.qMsgBox('Clusters already in order');
                else
                    jrclust.utils.qMsgBox(sprintf('%d clusters changed', nChanged));
                end
                obj.spatial_idx = val;
                [~,obj.channel_idx] = sort(obj.spatial_idx);                
                obj.updateFigWav();
                obj.updateSelect(1,1);                      
            end
        end

        % nShown
        function val = get.nShown(obj)
            val = numel(obj.showSubset);
        end

        % selected
        function set.selected(obj, val)
            if ~all(ismember(val, obj.showSubset))
                [~, argm] = min(abs(val' - obj.showSubset), [], 2);
                val = obj.showSubset(unique(argm));
            end

            obj.selected = val;
        end

        % showSubset
        function set.showSubset(obj, val)
            if isempty(val)
                val = 1:obj.hClust.nClusters;
            end

            obj.showSubset = val;
        end
    end
end
