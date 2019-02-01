classdef ControlFigure < dynamicprops & jrclust.views.Figure
    %CONTROLFIGURE
    
    properties
        hControls;
        hPanel;
    end

    %% LIFECYCLE
    methods
        function obj = ControlFigure(figTag, figPos, figName, figToolbar, figMenubar)
            obj = obj@jrclust.views.Figure(figTag, figPos, figName, figToolbar, figMenubar);

            obj.hPanel = uipanel('Parent', obj.hFig);
            obj.hControls = containers.Map();
        end
    end

    %% USER METHODS
    methods
        function axes(obj, varargin)
            %AXES Create new axes for this figure
            if ~obj.isReady
                obj.hFig = figure();
                obj.hPanel = uipanel('Parent', obj.hFig);
            end

            axes(obj.hPanel);
            obj.axApply('default', @hold, 'on');
        end

        function hCtl = addUicontrol(obj, controlKey, varargin)
            %ADDUICONTROL Create and store user interface control object
            if ~obj.isReady
                return;
            end

            obj.hControls(controlKey) = uicontrol('Parent', obj.hPanel, varargin{:});
            if nargout == 1
                hCtl = obj.hControls(controlKey);
            end
        end
    end

    %% GETTERS/SETTERS
    methods
    end
end

