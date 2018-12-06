classdef Figure < handle
    %FIGURE handle for JRCLUST manual figure
    %   Base class for specific figure types

    properties (Access=private, Hidden, SetObservable)
        hFig;           % Figure object
        hPlots;         % hashmap of current plots (return value of `plot`)
    end

    properties (Dependent, SetObservable)
        figTag;
        figPos;         
        figName;        % 
        figData;        % 'UserData' for hFig
        isReady;        % convenience prop for hFig validity
        outerPosition;  % position of the figure on the screen
    end

    properties (Hidden, SetObservable)
        hFunClick;      % left-click function handle
        hFunKey;        % keypress function handle
    end

    properties (Hidden, SetAccess=private)
        hideOnDrag;     % cell of plotKeys which we want to hide when we drag
        isMouseable;
        mouseStatus;
        prevPoint;
    end

    %% LIFECYCLE
    methods
        function obj = Figure(figTag, figPos, figName, figToolbar, figMenubar)
            %FIGURE Construct an instance of this class
            obj.hFig = figure();

            obj.figTag = figTag;
            obj.figPos = figPos;
            obj.figName = figName;

            set(obj.hFig, 'Color', 'w', 'NumberTitle', 'off');
            
            if figToolbar
                set(obj.hFig, 'ToolBar', 'figure');
            else
                set(obj.hFig, 'ToolBar', 'none');
            end

            if figMenubar
                set(obj.hFig, 'MenuBar', 'figure');
            else
                set(obj.hFig, 'MenuBar', 'none');
            end

            obj.hPlots = containers.Map();

            obj.hideOnDrag = {};
            obj.isMouseable = false;
        end
    end

    %% MOUSE METHODS
    methods (Access=private, Hidden)
        function scrollZoom(obj, varargin)
            %SCROLLZOOM Zoom with the scroll wheel
            hAx = get(obj.hFig, 'CurrentAxes');
            scrolls = varargin{2}.VerticalScrollCount;

            % get original limits
            xlim = get(hAx, 'xlim'); 
            ylim = get(hAx, 'ylim');

            % get the current camera position, and save the [z]-value
            camPosZ = get(hAx, 'CameraPosition'); 
            camPosZ = camPosZ(3);

            % get the current point
            oldPos = get(hAx, 'CurrentPoint');
            oldPos(1, 3) = camPosZ;

            % calculate zoom factor
            if scrolls < 0
                zoomFactor = sqrt(2); % 1 + 1/4;
            else
                zoomFactor = 1/sqrt(2); % 1 - 1/4;
            end

            set(hAx, 'CameraTarget', [oldPos(1, 1:2), 0], ...
                'CameraPosition', oldPos(1, 1:3));

            % adjust the camera view angle (equal to zooming in)
            camzoom(zoomFactor);
            keyFig = get(obj.hFig, 'CurrentCharacter');

            [zfX, zfY] = deal(zoomFactor); % zoomFactor X and Y
            if lower(keyFig) == 'x'
                zfY = 1;
            elseif lower(keyFig) == 'y'
                zfX = 1;
            end

            % zooming with the camera has the side-effect of
            % NOT adjusting the axes limits. We have to correct for this:
            xlim1 = (oldPos(1, 1) - min(xlim))/zfX;
            xlim2 = (max(xlim) - oldPos(1, 1))/zfX;
            xlim  = [oldPos(1,1) - xlim1, oldPos(1,1) + xlim2];
            set(hAx, 'xlim', xlim);

            ylim1 = (oldPos(1,2) - min(ylim))/zfY;
            ylim2 = (max(ylim) - oldPos(1,2))/zfY;
            ylim  = [oldPos(1,2) - ylim1, oldPos(1,2) + ylim2];            
            set(hAx, 'ylim', ylim);

            % set new camera position
            newPos = get(hAx, 'CurrentPoint');
            oldCamTarget =  get(hAx, 'CameraTarget');
            oldCamTarget(3) = camPosZ;
            newCamPos = oldCamTarget - newPos(1,1:3) + oldCamTarget(1,1:3);

            % adjust camera target and position
            set(hAx, 'CameraPosition', newCamPos(1, 1:3), ...
                'CameraTarget', [newCamPos(1, 1:2), 0]);

            % we also have to re-set the axes to stretch-to-fill mode
            set(hAx, 'CameraViewAngleMode', 'auto', ...
                'CameraPositionMode', 'auto', ...
                'CameraTargetMode', 'auto');
        end

        function panClick(obj, varargin)
            %PANCLICK Pan on mouse click
            hAx = get(obj.hFig, 'CurrentAxes');

            % select on click, pan on shift-click
            clickType = get(obj.hFig, 'SelectionType');

            if strcmp(clickType, 'normal') || strcmp(clickType, 'alt') % left or right click
                xyPoint = get(hAx, 'CurrentPoint');

                if isa(obj.hFunClick, 'function_handle')
                    obj.hFunClick(xyPoint([1, 3]), clickType);
                end
            elseif strcmp(clickType, 'extend') % shift-click
                obj.mouseStatus = 'down';
                obj.hideDrag(); % hide objects
                obj.prevPoint = get(hAx, 'CurrentPoint');
            end
        end

        function panMotion(obj, varargin)
            %PANMOUSE Move the mouse (with button clicked)
            if isempty(obj.prevPoint) || isempty(obj.mouseStatus)
                return
            end

            hAx = get(obj.hFig, 'CurrentAxes');
            % get current location (in pixels)
            curPoint = get(hAx, 'CurrentPoint');
            % get current XY-limits
            xlim = get(hAx, 'xlim');
            ylim = get(hAx, 'ylim');

            % find change in position
            dPoints = curPoint - obj.prevPoint;

            % adjust limits
            xlimNew = xlim - dPoints(1);
            ylimNew = ylim - dPoints(3);

            % set new limits
            set(hAx, 'xlim', xlimNew);
            set(hAx, 'ylim', ylimNew);

            % save new position
            obj.prevPoint = get(hAx, 'CurrentPoint');
        end

        function panRelease(obj, varargin)
            %PANRELEASE Release the mouse button after pan
            obj.showDrag();
            obj.mouseStatus = '';
        end

        function hideDrag(obj)
            %HIDEDRAG Hide objects on drag mouse
            if isempty(obj.hideOnDrag)
                return;
            end

            obj.toggleVisible(obj.hideOnDrag, 0);
        end

        function showDrag(obj)
            %SHOWDRAG Show figures after drag release
            try 
                obj.toggleVisible(obj.hideOnDrag, 1);
            catch
            end
        end

        function toggleVisible(obj, plotKey, fVis)
            %TOGGLEVISIBLE Toggle visibility of plot by key
            if isempty(plotKey)
                return;
            end

            if iscell(plotKey)
                if nargin < 3
                    cellfun(@(pKey) obj.toggleVisible(pKey), plotKey);
                else
                    cellfun(@(pKey) obj.toggleVisible(pKey, fVis), plotKey);
                end

                return;
            end

            if ~obj.hasPlot(plotKey)
                return;
            end

            hPlot = obj.hPlots(plotKey);
            if nargin == 2
                if strcmp(get(hPlot, 'Visible'), 'on')
                    set(hPlot, 'Visible', 'off');
                else
                    set(hPlot, 'Visible', 'on');
                end
            else
                if fVis == 0 % off
                    set(hPlot, 'Visible', 'off');
                elseif fVis == 1 % on
                    set(hPlot, 'Visible', 'on');
                end
            end

%             try
%                 if nargin == 2
%                     isVisible = false(size(plotKey));
%                     % toggle visibility
%                     for iPlot = 1:numel(plotKey)
%                         hPlot1 = plotKey(iPlot);
%                         if strcmpi(get(hPlot1, 'Visible'), 'on')
%                             isVisible(iPlot) = false;
%                             set(hPlot1, 'Visible', 'off');
%                         else
%                             isVisible(iPlot) = true;
%                             set(hPlot1, 'Visible', 'on');
%                         end
%                     end
%                 else
%                     % set visible directly
%                     if fVis
%                         isVisible = true(size(plotKey));
%                         set(plotKey, 'Visible', 'on');
%                     else
%                         isVisible = false(size(plotKey));
%                         set(plotKey, 'Visible', 'off');
%                     end
%                 end
%             catch
%                 return;
%             end
        end
    end

    %% USER METHODS
    methods
        function addDiag(obj, plotKey, lim, varargin)
            %ADDDIAG
            [xVals, yVals] = getDiagXY(lim);
            obj.addPlot(plotKey, xVals, yVals, varargin{:});
        end

        function addImagesc(obj, plotKey, varargin)
            %ADDLINE Create and store an image with scaled colors
            if obj.isReady
                hAx = obj.gca();
                if isempty(hAx)
                    obj.toForeground();
                    obj.hPlots(plotKey) = imagesc(varargin{:});
                else
                    obj.hPlots(plotKey) = imagesc(hAx, varargin{:});
                end

                if obj.isMouseable
                    obj.setMouseable();
                end
            end
        end

        function addLine(obj, plotKey, varargin)
            %ADDLINE Create and store a primitive line plot
            if obj.isReady
                hAx = obj.gca();
                if isempty(hAx)
                    obj.toForeground();
                    obj.hPlots(plotKey) = line(varargin{:});
                else
                    obj.hPlots(plotKey) = line(hAx, varargin{:});
                end

                if obj.isMouseable
                    obj.setMouseable();
                end
            end
        end

        function addPlot(obj, plotKey, varargin)
            %ADDPLOT Create and store a 2-D line plot
            if obj.isReady
                hAx = obj.gca();
                if isempty(hAx)
                    obj.toForeground();
                    obj.hPlots(plotKey) = plot(varargin{:});
                else
                    obj.hPlots(plotKey) = plot(hAx, varargin{:});
                end

                if obj.isMouseable
                    obj.setMouseable();
                end
            end
        end

        function axes(obj, varargin)
            %AXES Create new axes for this figure
            if ~obj.isReady
                obj.hFig = figure();
            end

            obj.clf();
            axes(obj.hFig);
            obj.hold('on');
        end

        function axis(obj, varargin)
            %AXIS Set axis limits and aspect ratios
            hAx = obj.gca();
            if ~isempty(hAx)
                axis(hAx, varargin{:});
            end
        end

        function val = axGet(obj, varargin)
            %AXGET Query current axes properties
            hAx = obj.gca();
            if ~isempty(hAx)
                val = get(hAx, varargin{:});
            else
                val = [];
            end
        end

        function axSet(obj, varargin)
            %AXSET Set current axes properties
            hAx = obj.gca();
            if ~isempty(hAx)
                set(hAx, varargin{:});
            end
        end

        function clf(obj)
            %CLF Clear figure window and remove references to plots
            if obj.isReady
                clf(obj.hFig);
                obj.hPlots = containers.Map(); % reset mapping
            end
        end

        function close(obj)
            %CLOSE Close the figure
            if obj.isReady
                close(obj.hFig);
            end
        end

        function colorbar(obj, varargin)
            %COLORBAR Colorbar showing color scale
            if obj.isReady
                hAx = get(obj.hFig, 'CurrentAxes');
                colorbar(hAx);
            end
        end

        function val = figGet(obj, varargin)
            %FIGGET Query figure properties
            if obj.isReady
                val = get(obj.hFig, varargin{:});
            else
                val = [];
            end
        end

        function figSet(obj, varargin)
            %FIGSET Set figure properties
            if obj.isReady
                set(obj.hFig, varargin{:});
            end
        end

        function grid(obj, varargin)
            if obj.isReady
                hAx = get(obj.hFig, 'CurrentAxes');
                grid(hAx, varargin{:});
            end
        end

        function hp = hasPlot(obj, plotKey)
            %HASPLOT Return true iff a plot exists with plotKey as label
            hp = ischar(plotKey) && isKey(obj.hPlots, plotKey);
        end

        function hidePlot(obj, plotKey)
            %HIDEPLOT Set XData and YData of a plot to nan
            obj.updatePlot(plotKey, nan, nan);
        end

        function hold(obj, varargin)
            %HOLD Retain current plot when adding new plots
            if obj.isReady
                hAx = get(obj.hFig, 'CurrentAxes');
                hold(hAx, varargin{:});
            end
        end

        function rmPlot(obj, plotKey)
            %RMPLOT Remove a plot by key
            if ~obj.hasPlot(plotKey)
                return;
            end

            delete(obj.hPlots(plotKey));
            remove(obj.hPlots, plotKey);
        end

        function setHideOnDrag(obj, plotKey)
            if ~obj.hasPlot(plotKey) || any(strcmp(obj.hideOnDrag, plotKey))
                return;
            end

            obj.hideOnDrag{end+1} = plotKey;
        end

        function setMouseable(obj, hFunClick)
            %SETMOUSEABLE Set the figure to be mouseable
            obj.isMouseable = true;
            if nargin == 2 && ~isempty(hFunClick)
                obj.hFunClick = hFunClick;
            end                

            if obj.isReady
                hAx = get(obj.hFig, 'CurrentAxes');
                 % is2D might disappear in a future release...
                if ~is2D(hAx) || isempty(hAx)
                    return;
                end

                % define zoom with scroll wheel, pan with left click
                set(obj.hFig, ...
                    'WindowScrollWheelFcn' , @obj.scrollZoom, ...
                    'WindowButtonDownFcn'  , @obj.panClick, ...
                    'WindowButtonUpFcn'    , @obj.panRelease, ...
                    'WindowButtonMotionFcn', @obj.panMotion);
            end
        end

        function title(obj, t)
            if obj.isReady
                hAx = get(obj.hFig, 'CurrentAxes');
                title(hAx, t, 'Interpreter', 'none', 'FontWeight', 'normal');
            end
        end

        function toForeground(obj)
            %TOFOREGROUND Move current plot into focus
            if obj.isReady
                figure(obj.hFig);
            end
        end

        function hMenu = uimenu(obj, varargin)
            if obj.isReady
                hMenu = uimenu(obj.hFig, varargin{:});
            else
                hMenu = [];
            end
        end

        function updateImagesc(obj, plotKey, CData)
            %UPDATEIMAGESC Set CData of an image
            if ~obj.hasPlot(plotKey)
                return;
            end

            hPlot = obj.hPlots(plotKey);
            set(hPlot, 'CData', CData);
        end

        function updatePlot(obj, plotKey, newXData, newYData, UserData)
            %UPDATEPLOT Set XData, YData, and optionally UserData of a plot
            if ~obj.hasPlot(plotKey)
                return;
            end

            % clear data if we're empty
            if isempty(newXData) || isempty(newYData)
                obj.hidePlot(plotKey);
                return;
            end

            if nargin < 4
                UserData = [];
            end

            hPlot = obj.hPlots(plotKey);

            % only update if both x and y are changed
            oldXData = get(hPlot, 'XData');
            oldYData = get(hPlot, 'YData');

            doUpdate = true;
            if (numel(oldXData) == numel(newXData)) && (numel(oldYData) == numel(newYData))
                if (std(oldXData(:) - newXData(:)) == 0) && (std(oldYData(:) - newYData(:)) == 0)
                    doUpdate = false;
                end
            end

            if doUpdate
                set(hPlot, 'XData', newXData, 'YData', newYData);
            end
            if ~isempty(UserData)
                set(hPlot, 'UserData', UserData);
            end
        end

        function xlabel(obj, varargin)
            if obj.isReady
                hAx = get(obj.hFig, 'CurrentAxes');
                xlabel(hAx, varargin{:});
            end
        end

        function ylabel(obj, varargin)
            if obj.isReady
                hAx = get(obj.hFig, 'CurrentAxes');
                ylabel(hAx, varargin{:});
            end
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function hAx = gca(obj)
            %GCA Get current axes
            if obj.isReady
                hAx = get(obj.hFig, 'CurrentAxes');
            else
                hAx = [];
            end
        end
    end

    %% GETTERS/SETTERS
    methods
        % figData
        function fd = get.figData(obj)
            if obj.isReady
                fd = get(obj.hFig, 'UserData');
            else
                fd = [];
            end
        end
        function set.figData(obj, fm)
            if obj.isReady
                set(obj.hFig, 'UserData', fm);
            end
        end

        % figName
        function fn = get.figName(obj)
            if obj.isReady
                fn = get(obj.hFig, 'Name');
            else
                fn = '';
            end
        end
        function set.figName(obj, figName)
            if obj.isReady
                set(obj.hFig, 'Name', figName);
            end
        end

        % figPos
        function fp = get.figPos(obj)
            if obj.isReady
                fp = get(obj.hFig, 'OuterPosition');
            else
                fp = [];
            end
        end
        function set.figPos(obj, figPos)
            if isempty(figPos)
                return;
            end

            tbHeight = 40; % taskbar height

            rootPos = get(groot, 'ScreenSize');
            screenWidth = rootPos(3);
            screenHeight = rootPos(4) - tbHeight;

            figPos = [max(1, round(figPos(1) * screenWidth)) ... % left
                      tbHeight + max(1, round(figPos(2)*screenHeight)) ... % bottom
                      min(screenWidth, round(figPos(3)*screenWidth)) ... % right
                      min(screenHeight, round(figPos(4)*screenHeight))]; % top

            drawnow;
            set(obj.hFig, 'OuterPosition', figPos);
        end

        % figTag
        function ft = get.figTag(obj)
            if obj.isReady
                ft = get(obj.hFig, 'Tag');
            else
                ft = '';
            end
        end
        function set.figTag(obj, figTag)
            if obj.isReady
                set(obj.hFig, 'Tag', figTag);
            end
        end

        % hFunClick
        function set.hFunClick(obj, hf)
            failMsg = 'hFunClick must be a function handle';
            assert(isa(hf, 'function_handle'), failMsg);
            obj.hFunClick = hf;
        end

        % hFunKey
        function set.hFunKey(obj, hf)
            failMsg = 'hFunKey must be a function handle';
            assert(isa(hf, 'function_handle'), failMsg);
            obj.hFunKey = hf;
            if obj.isReady
                set(obj.hFig, 'KeyPressFcn', obj.hFunKey);
            end
        end

        % isReady
        function ir = get.isReady(obj)
            ir = (~isempty(obj.hFig) && isvalid(obj.hFig));
        end

        % outerPosition
        function op = get.outerPosition(obj)
            if obj.isReady
                op = get(obj.hFig, 'OuterPosition');
            else
                op = [];
            end
        end
        function set.outerPosition(obj, op)
            if obj.isReady
                set(obj.hFig, 'OuterPosition', op);
            end
        end
    end
end

