classdef CurateController < handle
    %CURATECONTROLLER Interface for manually curating sorted clusters
    properties (SetAccess=private, Hidden, SetObservable)
        cRes;           % curate results struct, returned at endSession
        hClust;         % Clustering object
    end

    properties (Dependent)
        hCfg;
    end

    properties (AbortSet, SetAccess=private, Hidden, Transient, SetObservable)
        currentSite;    % current site in FigTime
        hFigs;          % containers.Map of Figure objects
        hMenus;         % containers.Map of Menu handles
        isEnding;       % to prevent double prompting when endSession is called
        isWorking;      % to prevent keypress/mouseclick functions from colliding
        maxAmp;         % scaling factor for 
        projSites;      % current sites in FigProj
        selected;       % selected clusters, in order of selection
    end

    %% LIFECYCLE
    methods
        function obj = CurateController(hClust)
            obj.hClust = hClust;
            obj.hFigs = containers.Map();
            obj.hMenus = containers.Map();
            obj.isEnding = 0;
            obj.isWorking = 0;
            obj.currentSite = [];
            obj.maxAmp = [];
            obj.selected = [];
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
    methods
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
        annotateUnit(obj, note, doConfirm);
        autoDelete(obj);
        autoMerge(obj);
        autoSplit(obj, multisite);
        closeFigures(obj);
        deleteAnnotated(obj);
        deleteClusters(obj, deleteMe);
        exportFiringRate(obj);
        exportMeanWf(obj, exportAll);
        exportTraces(obj);
        res = figApply(obj, hFun, varargin);
        killFigWav(obj, hObject, hEvent);
        mergeSelected(obj);
        plotAllFigures(obj);
        spawnFigures(obj);
        splitCluster(obj, iCluster, retainedSpikes);
        toggleRaw(obj, hMenu);
        updateMenu(obj);
        updateProjection(obj, proj);
        updateSelect(obj, iClusters);
    end
    

    %% GETTERS/SETTERS
    methods
        % hCfg
        function hc = get.hCfg(obj)
            if ~isempty(obj.hClust)
                hc = obj.hClust.hCfg;
            else
                hc = [];
            end
        end

        % hClust
        function hc = get.hClust(obj)
            if ~isempty(obj.cRes) && isfield(obj.cRes, 'hClust')
                hc = obj.cRes.hClust;
            else
                hc = obj.hClust;
            end
        end
    end
end
