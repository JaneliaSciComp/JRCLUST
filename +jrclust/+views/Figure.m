classdef Figure < handle
    %FIGURE handle for JRCLUST manual figure
    %   Base class for specific figure types

    properties (Access=protected, Hidden, SetObservable)
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

    properties (Hidden, SetAccess=protected)
        hideOnDrag;     % cell of plotKeys which we want to hide when we drag
        isMouseable;
        isWaiting;      % are we waiting on an operation to complete?
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
            obj.isWaiting = false;
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

    end

    %% USER METHODS
    methods
        function addBar(obj, plotKey, varargin)
            %ADDBAR Create and store a bar graph
            if obj.isReady
                hAx = obj.gca();
                if isempty(hAx)
                    obj.toForeground();
                    obj.hPlots(plotKey) = bar(varargin{:});
                else
                    obj.hPlots(plotKey) = bar(hAx, varargin{:});
                end

                if obj.isMouseable
                    obj.setMouseable();
                end
            end
        end

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

        function addPatch(obj, plotKey, varargin)
            %ADDPATCH Create and store one or more filled polygons
            if obj.isReady
                hAx = obj.gca();
                if isempty(hAx)
                    obj.toForeground();
                    obj.hPlots(plotKey) = patch(varargin{:});
                else
                    obj.hPlots(plotKey) = patch(hAx, varargin{:});
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

        function addStairs(obj, plotKey, varargin)
            %ADDSTAIRS Create and store a stairstep graph
            if obj.isReady
                hAx = obj.gca();
                if isempty(hAx)
                    obj.toForeground();
                    obj.hPlots(plotKey) = stairs(varargin{:});
                else
                    obj.hPlots(plotKey) = stairs(hAx, varargin{:});
                end

                if obj.isMouseable
                    obj.setMouseable();
                end
            end
        end

        function addText(obj, plotKey, varargin)
            %ADDSTAIRS Add and store text descriptions for data points
            if obj.isReady
                hAx = obj.gca();
                if isempty(hAx)
                    obj.toForeground();
                    obj.hPlots(plotKey) = text(varargin{:});
                else
                    obj.hPlots(plotKey) = text(hAx, varargin{:});
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

        function caxis(obj, limits)
            %CAXIS Set colormap limits for current axis
            hAx = obj.gca();
            if ~isempty(hAx)
                caxis(hAx, limits);
            end
        end

        function cla(obj)
            %CLA Clear axes
            hAx = obj.gca();
            if ~isempty(hAx)
                cla(hAx);
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
                hAx = obj.gca();
                hold(hAx, varargin{:});
            end
        end

        function [plotKey, yOffsets] = multiplot(obj, plotKey, scale, XData, YData, yOffsets, fScatter)
            % Create (nargin > 2) or rescale (nargin <= 2) a multi-line plot
            % TODO: separate these (presumably multiple plots are passed in somewhere)
            if nargin <= 2 % rescale
                obj.rescalePlot(plotKey, scale);
                % handle_fun_(@rescale_plot_, plotKey, scale);
                yOffsets = [];
                return;
            end

            if nargin < 7
                fScatter = false;
            end

            if isa(YData, 'gpuArray')
                YData = gather(YData);
            end
            if isa(YData, 'int16')
                YData = single(YData);
            end

            if nargin < 6
                yOffsets = 1:size(YData, 2);
            end

            [plotKey, yOffsets] = doMultiplot(obj, plotKey, scale, XData, YData, yOffsets, fScatter);
        end

        function val = plotGet(obj, plotKey, varargin)
            %PLOTGET Query a plot by key
            if ~obj.hasPlot(plotKey)
                return;
            end

            val = get(obj.hPlots(plotKey), varargin{:});
        end

        function plotSet(obj, plotKey, varargin)
            %PLOTSET Set a plot value by key
            if ~obj.hasPlot(plotKey)
                return;
            end

            set(obj.hPlots(plotKey), varargin{:});
        end

        function rescalePlot(obj, plotKey, scale)
            % hPlot must have UserData containing scale, shape
            if ~obj.hasPlot(plotKey)
                return;
            end

            hPlot = obj.hPlots(plotKey);
            userData = get(hPlot, 'UserData');

            % backward compatible
            % if ~isfield(userData, 'yOffsets')
            %     userData.yOffsets = 1:userData.shape(2);
            % end

            % Rescale
            YData = reshape(get(hPlot, 'YData'), userData.shape);
            if isfield(userData, 'fScatter')
                fScatter = userData.fScatter; 
            else
                fScatter = false;
            end

            if fScatter
                YData = (YData(:) - userData.yOffsets(:)) * userData.scale; %restore original 
                YData = YData(:) / scale + userData.yOffsets(:); %convert to plot
            elseif isvector(userData.yOffsets)
                YData = bsxfun(@minus, YData, userData.yOffsets(:)') * userData.scale; % restore original 
                YData = bsxfun(@plus, YData/scale, userData.yOffsets(:)'); % convert to plot
            else
                for iSpike = 1:userData.shape(3)
                    iYOffsets = userData.yOffsets(:, iSpike)';
                    iYData = bsxfun(@minus, YData(:, :, iSpike), iYOffsets)*userData.scale;
                    YData(:, :, iSpike) = bsxfun(@plus, iYData/scale, iYOffsets);
                end  
            end

            userData.scale = scale; % update local scale
            set(hPlot, 'YData', YData(:), 'UserData', userData);
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

        function setWindow(obj, xlim1, ylim1, xlim0, ylim0)
            %SETWINDOW set the window within the box limit
            if nargin <= 3
                % square case
                xlim0 = ylim1;
                ylim1 = xlim1;
                ylim0 = xlim0;
            end

            lastFocused = gcf();
            obj.toForeground();

            dx = diff(xlim1);
            dy = diff(ylim1);

%             if xlim1(1)<xlim0(1), xlim1 = xlim0(1) + [0, dx]; end
%             if xlim1(2)>xlim0(2), xlim1 = xlim0(2) + [-dx, 0]; end
%             if ylim1(1)<ylim0(1), ylim1 = ylim0(1) + [0, dy]; end
%             if ylim1(2)>ylim0(2), ylim1 = ylim0(2) + [-dy, 0]; end
% 
%             xlim1(1) = max(xlim1(1), xlim0(1));
%             ylim1(1) = max(ylim1(1), ylim0(1));
%             xlim1(2) = min(xlim1(2), xlim0(2));
%             ylim1(2) = min(ylim1(2), ylim0(2));

            xlim1 = jrclust.utils.trimLim(xlim1, xlim0);
            ylim1 = jrclust.utils.trimLim(ylim1, ylim0);

            obj.axis([xlim1, ylim1]);

            figure(lastFocused);
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

        function uistack(obj, plotKey, varargin)
            %UISTACK Reorder visual stacking order of plot given by plotKey
            if ~obj.hasPlot(plotKey)
                return;
            end

            hPlot = obj.hPlots(plotKey);
            uistack(hPlot, varargin{:});
        end

        function updateImagesc(obj, plotKey, CData)
            %UPDATEIMAGESC Set CData of an image
            if ~obj.hasPlot(plotKey)
                return;
            end

            hPlot = obj.hPlots(plotKey);
            set(hPlot, 'CData', CData);
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

            if nargin < 5
                UserData = [];
            end

            hPlot = obj.hPlots(plotKey);

            % only update if both x and y are changed
            oldXData = get(hPlot, 'XData');
            oldYData = get(hPlot, 'YData');

            doUpdate = true;
            if (numel(oldXData) == numel(newXData)) && (numel(oldYData) == numel(newYData))
                if (norm(oldXData(:) - newXData(:)) == 0) && (norm(oldYData(:) - newYData(:)) == 0)
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

        function wait(obj, doWait)
            %WAIT Toggle waiting status and set pointer to reflect
            if obj.isReady
                if nargin == 2 && doWait
                    obj.isWaiting = true;
                elseif nargin == 2 % ~doWait
                    obj.isWaiting = false;
                else % nargin == 0, toggle on/off
                    obj.isWaiting = ~obj.isWaiting;
                end

                if obj.isWaiting
                    obj.figSet('Pointer', 'watch');
                else
                    obj.figSet('Pointer', 'arrow');
                end
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
                if isempty(hAx)
                    obj.axes();
                    hAx = get(obj.hFig, 'CurrentAxes');
                end
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

