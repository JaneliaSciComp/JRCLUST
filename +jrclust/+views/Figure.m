classdef Figure < handle
    %FIGURE handle for JRCLUST manual figure
    %   Base class for specific figure types
    properties (SetAccess=protected, Hidden, SetObservable)
        hFig;           % Figure object
        hAxes;          % hashmap of current axes (return value of `axes`)
        hPlots;         % hashmap of current plots (return value of `plot`)
        hSubplots;      % hashmap of current subplots (return value of `subplot`)
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
        initialPos;
    end

    properties (Hidden, SetAccess=protected)
        hideOnDrag;     % plotKeys to hide when we drag
        hideOnDragOff;  % plotKeys normally hidden on drag and manually toggled off
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
            obj.initialPos = figPos;

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

            obj.hAxes = containers.Map();
            obj.hPlots = containers.Map();
            obj.hSubplots = containers.Map();

            obj.hideOnDrag = {};
            obj.hideOnDragOff = {};
            obj.isMouseable = 0;
            obj.isWaiting = 0;
        end
    end

    %% MOUSE METHODS
    methods (Access=private, Hidden)
        function scrollZoom(obj, axKey, hEvent)
            %SCROLLZOOM Zoom with the scroll wheel
            hAx = obj.gca(axKey);
            scrolls = hEvent.VerticalScrollCount;

            % get original limits
            XLim = get(hAx, 'XLim'); 
            YLim = get(hAx, 'YLim');

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
            xlim1 = (oldPos(1, 1) - min(XLim))/zfX;
            xlim2 = (max(XLim) - oldPos(1, 1))/zfX;
            XLim  = [oldPos(1,1) - xlim1, oldPos(1,1) + xlim2];
            set(hAx, 'XLim', XLim);

            ylim1 = (oldPos(1,2) - min(YLim))/zfY;
            ylim2 = (max(YLim) - oldPos(1,2))/zfY;
            YLim  = [oldPos(1,2) - ylim1, oldPos(1,2) + ylim2];            
            set(hAx, 'YLim', YLim);

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

        function panClick(obj, axKey)
            %PANCLICK Pan on mouse click
            hAx = obj.gca(axKey);

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

        function panMotion(obj, axKey)
            %PANMOUSE Move the mouse (with button clicked)
            if isempty(obj.prevPoint) || isempty(obj.mouseStatus)
                return
            end

            hAx = obj.gca(axKey);
            % get current location (in pixels)
            curPoint = get(hAx, 'CurrentPoint');
            % get current XY-limits
            XLim = get(hAx, 'XLim');
            YLim = get(hAx, 'YLim');

            % find change in position
            dPoints = curPoint - obj.prevPoint;

            % adjust limits
            XlimNew = XLim - dPoints(1);
            YLimNew = YLim - dPoints(3);

            % set new limits
            set(hAx, 'XLim', XlimNew);
            set(hAx, 'YLim', YLimNew);

            % save new position
            obj.prevPoint = get(hAx, 'CurrentPoint');
        end

        function panRelease(obj)
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
        function hAx = addAxes(obj, axKey, varargin)
            %ADDAXES Create and store new axes for this figure
            if ~obj.isReady
                obj.hFig = figure();
            end

            if obj.hasAxes(axKey)
                return;
            end

            hAx = axes(obj.hFig, varargin{:});
            obj.hAxes(axKey) = hAx;
            obj.axApply(axKey, @hold, 'on');
        end

        function addDiag(obj, plotKey, lim, varargin)
            %ADDDIAG
            [XVals, YVals] = getDiagXY(lim);
            obj.addPlot(plotKey, XVals, YVals, varargin{:});
        end

        function addPlot(obj, plotKey, hFunPlot, varargin)
            %ADDPLOT Create and store a plot
            if ~isa(hFunPlot, 'function_handle')
                varargin = [hFunPlot, varargin];
                hFunPlot = @plot;
            end

            if obj.isReady
                if ~isempty(varargin) && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                    hAx = varargin{1}; % varargin{1} is an Axes
                    varargin = varargin(2:end);
                else % add to default axis
                    hAx = obj.gca();
                end

                if obj.hasPlot(plotKey)
                    obj.rmPlot(plotKey);
                end
                obj.hPlots(plotKey) = hFunPlot(hAx, varargin{:});
            end
        end

        function addSubplot(obj, plotKey, nRows, nCols)
            %ADDSUBPLOTGRID Create and store axes in tiled positions
            if obj.isReady
                sp = matlab.graphics.axis.Axes.empty(nRows*nCols, 0);
                for i = 1:nRows*nCols
                    sp(i) = subplot(nRows, nCols, i);
                end
                obj.hSubplots(plotKey) = sp;
            end
        end

        function addTable(obj, plotKey, lim, varargin)
            %ADDTABLE Add a table
            XData = floor((lim(1)*2:lim(2)*2+1)/2);
            YData = repmat([lim(1), lim(2), lim(2), lim(1)], [1, ceil(numel(XData)/4)]);
            YData = YData(1:numel(XData));
            obj.addPlot(plotKey, [XData(1:end-1), fliplr(YData)], [YData(1:end-1), fliplr(XData)], varargin{:});
        end

        function vals = axApply(obj, axKey, hFun, varargin)
            %AXAPPLY Apply a function to current axes
            if ~obj.hasAxes(axKey)
                return;
            end

            hAx = obj.hAxes(axKey);
            if ~isempty(hAx) && isa(hFun, 'function_handle')
                if strcmp(func2str(hFun), 'title') % default arguments to title
                    varargin = [varargin, {'FontWeight', 'normal'}];
                end

                % apply hFun
                if nargout == 1
                    vals = hFun(hAx, varargin{:});
                else
                    hFun(hAx, varargin{:});
                end
            else
                vals = [];
            end
        end

        function cla(obj)
            %CLA Clear axes and remove all plots
            cellfun(@(axKey) cla(obj.hAxes(axKey)), keys(obj.hAxes));
            obj.hAxes = containers.Map();
            obj.hPlots = containers.Map();
            obj.hSubplots = containers.Map();
        end

        function clf(obj)
            %CLF Clear figure window and remove references to plots
            if obj.isReady
                clf(obj.hFig);

                 % reset mappings
                obj.hAxes = containers.Map();
                obj.hPlots = containers.Map();
                obj.hSubplots = containers.Map();
            end
        end

        function close(obj)
            %CLOSE Close the figure
            if obj.isReady
                jrclust.utils.tryClose(obj.hFig);
            end
        end

        function vals = figApply(obj, hFun, varargin)
            %FIGAPPLY Apply a function to hFig
            if obj.isReady && isa(hFun, 'function_handle')
                if nargout > 0
                    vals = hFun(obj.hFig, varargin{:});
                else
                    hFun(obj.hFig, varargin{:});
                end
            else
                vals = [];
            end
        end

        function hp = hasAxes(obj, axKey)
            %HASAXES Return true iff a plot exists with plotKey as label
            hp = ischar(axKey) && isKey(obj.hAxes, axKey);
        end

        function hp = hasPlot(obj, plotKey)
            %HASPLOT Return true iff a plot exists with plotKey as label
            hp = ischar(plotKey) && isKey(obj.hPlots, plotKey);
        end

        function hp = hasSubplot(obj, plotKey)
            %HASSUBPLOT Return true iff a subplot exists with plotKey as label
            hp = ischar(plotKey) && isKey(obj.hSubplots, plotKey);
        end

        function hidePlot(obj, plotKey)
            %HIDEPLOT Set XData and YData of a plot to nan
            obj.updatePlot(plotKey, nan, nan); % updatePlot checks for existence of plotKey
        end

        function [plotKey, YOffsets] = multiplot(obj, plotKey, scale, XData, YData, YOffsets, fScatter)
            % Create (nargin > 2) or rescale (nargin <= 2) a multi-line plot
            % TODO: separate these (presumably multiple plots are passed in somewhere)
            if nargin <= 2 % rescale
                obj.rescalePlot(plotKey, scale);
                % handle_fun_(@rescale_plot_, plotKey, scale);
                YOffsets = [];
                return;
            end

            if nargin < 7
                fScatter = 0;
            end

            if isa(YData, 'gpuArray')
                YData = gather(YData);
            end
            if isa(YData, 'int16')
                YData = single(YData);
            end

            if nargin < 6
                YOffsets = 1:size(YData, 2);
            end

            [plotKey, YOffsets] = doMultiplot(obj, plotKey, scale, XData, YData, YOffsets, fScatter);
        end

        function vals = plotApply(obj, plotKey, hFun, varargin)
            %PLOTAPPLY Apply hFun to the plot given by plotKey
            if ~obj.hasPlot(plotKey)
                return;
            end

            try
                if nargout == 1
                    vals = hFun(obj.hPlots(plotKey), varargin{:});
                else
                    hFun(obj.hPlots(plotKey), varargin{:});
                end
            catch
            end
        end

        function rescalePlot(obj, plotKey, scale)
            % hPlot must have UserData containing scale, shape
            if ~obj.hasPlot(plotKey)
                return;
            end

            hPlot = obj.hPlots(plotKey);
            UserData = get(hPlot, 'UserData');
            if ~isfield(UserData, 'scale') || ~isfield(UserData, 'shape')
                return;
            end

            % rescale
            YData = reshape(get(hPlot, 'YData'), UserData.shape);
            if isfield(UserData, 'fScatter')
                fScatter = UserData.fScatter; 
            else
                fScatter = 0;
            end

            if fScatter
                YData = (YData(:) - UserData.yOffsets(:)) * UserData.scale; %restore original 
                YData = YData(:) / scale + UserData.yOffsets(:); %convert to plot
            elseif isvector(UserData.yOffsets)
                YData = bsxfun(@minus, YData, UserData.yOffsets(:)') * UserData.scale; % restore original 
                YData = bsxfun(@plus, YData/scale, UserData.yOffsets(:)'); % convert to plot
            else
                for iSpike = 1:UserData.shape(3)
                    iYOffsets = UserData.yOffsets(:, iSpike)';
                    iYData = bsxfun(@minus, YData(:, :, iSpike), iYOffsets)*UserData.scale;
                    YData(:, :, iSpike) = bsxfun(@plus, iYData/scale, iYOffsets);
                end  
            end

            UserData.scale = scale; % update local scale
            set(hPlot, 'YData', YData(:), 'UserData', UserData);
        end

        function resetPos(obj)
            %RESETPOS Set figure position to where it was spawned initially
            obj.figPos = obj.initialPos;
        end

        function rmPlot(obj, plotKey)
            %RMPLOT Remove a plot by key
            if ~obj.hasPlot(plotKey)
                return;
            end

            delete(obj.hPlots(plotKey)); % delete the plot
            remove(obj.hPlots, plotKey); % remove the key
        end

        function setHideOnDrag(obj, plotKey)
            if ~obj.hasPlot(plotKey) || any(strcmp(obj.hideOnDrag, plotKey))
                return;
            end

            obj.hideOnDrag{end+1} = plotKey;
        end

        function setMouseable(obj, hFunClick, axKey)
            %SETMOUSEABLE Set the figure to be mouseable
            if nargin >= 2 && ~isempty(hFunClick)
                obj.hFunClick = hFunClick;
            end

            if nargin < 3
                axKey = 'default';
            end

            obj.isMouseable = 1;
            if obj.isReady && obj.hasAxes(axKey)
                % is2D might disappear in a future release...
                hAx = obj.hAxes(axKey);
                if ~is2D(hAx) || isempty(hAx)
                    return;
                end

                % define zoom with scroll wheel, pan with left click
                set(obj.hFig, ...
                    'WindowScrollWheelFcn' , @(hO, hE) obj.scrollZoom(axKey, hE), ...
                    'WindowButtonDownFcn'  , @(hO, hE) obj.panClick(axKey), ...
                    'WindowButtonUpFcn'    , @(hO, hE) obj.panRelease(), ...
                    'WindowButtonMotionFcn', @(hO, hE) obj.panMotion(axKey));
            end
        end

        function setWindow(obj, xlim1, ylim1, xlim0, ylim0)
            %SETWINDOW set the window within the box limit
            if nargin <= 3 % square case
                xlim0 = ylim1;
                ylim1 = xlim1;
                ylim0 = xlim0;
            end

            lastFocused = gcf();
            obj.toForeground();

            xlim1 = jrclust.utils.trimLim(xlim1, xlim0);
            ylim1 = jrclust.utils.trimLim(ylim1, ylim0);

            obj.axApply('default', @axis, [xlim1, ylim1]);

            figure(lastFocused);
        end

        function vals = subplotApply(obj, plotKey, spIndex, hFun, varargin)
            %SUBPLOTAPPLY Apply a plotting function to a subplot
            vals = [];

            if ~obj.hasSubplot(plotKey) || numel(obj.hSubplots(plotKey)) < spIndex
                return;
            end

            hAxes = obj.hSubplots(plotKey);
            hAx = hAxes(spIndex);

            % apply hFun
            if nargout == 1
                vals = hFun(hAx, varargin{:});
            else
                hFun(hAx, varargin{:});
            end

            % save hAxes
            obj.hSubplots(plotKey) = hAxes;
        end

        function toForeground(obj)
            %TOFOREGROUND Move current plot into focus
            if obj.isReady
                figure(obj.hFig);
            end
        end

        function fVis = toggleVisible(obj, plotKey, fVis)
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
                    fVis = 0;
                else
                    set(hPlot, 'Visible', 'on');
                    fVis = 1;
                end
            else
                if fVis == 0 % off
                    set(hPlot, 'Visible', 'off');
                elseif fVis == 1 % on
                    set(hPlot, 'Visible', 'on');
                end
            end
        end

        function hMenu = uimenu(obj, varargin)
            if obj.isReady
                hMenu = uimenu(obj.hFig, varargin{:});
            else
                hMenu = [];
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

            doUpdate = 1;
            if (numel(oldXData) == numel(newXData)) && (numel(oldYData) == numel(newYData))
                if all(oldXData(:) - newXData(:) == 0) && all(oldYData(:) - newYData(:) == 0)
                    doUpdate = 0;
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
                    obj.isWaiting = 1;
                elseif nargin == 2 % ~doWait
                    obj.isWaiting = 0;
                else % nargin == 0, toggle on/off
                    obj.isWaiting = ~obj.isWaiting;
                end

                if obj.isWaiting
                    obj.figApply(@set, 'Pointer', 'watch');
                else
                    obj.figApply(@set, 'Pointer', 'arrow');
                end
            end
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function hAx = gca(obj, axKey)
            %GCA Get current axes, optionally specifying a key
            if nargin < 2
                axKey = 'default';
            end

            if obj.isReady 
                if ~obj.hasAxes(axKey)
                    hAx = obj.addAxes(axKey);
                else
                    hAx = obj.hAxes(axKey);
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
            if obj.isReady %#ok<*MCSUP>
                set(obj.hFig, 'KeyPressFcn', obj.hFunKey);
            end
        end

        % isReady
        function ir = get.isReady(obj)
            ir = (~isempty(obj.hFig) && ishandle(obj.hFig));
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

