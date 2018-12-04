classdef Figure < handle
    %FIGURE handle for JRCLUST manual figure
    %   Base class for specific figure types

    properties (Access=private, Hidden, SetObservable)
        hFig;
    end

    properties (Dependent, SetObservable)
        figTag;
        figPos;         
        figName;        % 
        figData;        % 
        figMetadata;    % 'UserData' for hFig
        outerPosition;  % position of the figure on the screen
    end

    properties (Hidden, SetObservable)
        hFunClick;      % left-click function handle
        hFunKey;        % keypress function handle
        status;         % 
    end

    properties (Hidden, SetAccess=private)
        isMouseable;
        mouseFigHidden;
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

            obj.isMouseable = false;
            obj.mouseFigHidden = false;
            obj.status = 'ready';
            set(obj.hFig, 'CloseRequestFcn', @obj.destroyFig);
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
            if ~isfield(obj.figMetadata, 'cvhHide_mouse')
                return;
            end

            try
                vhHide = obj.figMetadata.cvhHide_mouse;
                if strcmpi(get(vhHide{1}(1), 'Visible'), 'off')
                    return;
                end

                obj.toggleVisible(vhHide, 0); % hide
                obj.mouseFigHidden = vhHide; % hidden object
            catch
            end
        end
        
        function showDrag(obj)
            if isempty(obj.mouseFigHidden)
                return;
            end 
            try 
                obj.toggleVisible(obj.mouseFigHidden, 1); % hide
                obj.mouseFigHidden = [];
            catch
            end
        end
        
        function vlVisible = toggleVisible(obj, vhPlot, fVisible)
            if isempty(vhPlot)
                return;
            end

            if iscell(vhPlot)
                cvhPlot = vhPlot;

                if nargin < 3
                    cellfun(@(vhPlot) obj.toggleVisible(vhPlot), cvhPlot);
                else
                    cellfun(@(vhPlot) obj.toggleVisible(vhPlot, fVisible), cvhPlot);
                end

                return;
            end

            try
                if nargin == 2
                    vlVisible = false(size(vhPlot));
                    % toggle visibility
                    for iH = 1:numel(vhPlot)
                        hPlot1 = vhPlot(iH);
                        if strcmpi(get(hPlot1, 'Visible'), 'on')
                            vlVisible(iH) = 0;
                            set(hPlot1, 'Visible', 'off');
                        else
                            vlVisible(iH) = 1;
                            set(hPlot1, 'Visible', 'on');
                        end
                    end
                else
                    % set visible directly
                    if fVisible
                        vlVisible = true(size(vhPlot));
                        set(vhPlot, 'Visible', 'on');
                    else
                        vlVisible = false(size(vhPlot));
                        set(vhPlot, 'Visible', 'off');
                    end
                end
            catch
                return;
            end
        end
    end
    
    %% USER METHODS
    methods
        function axes(obj, varargin)
            %AXES Create new axes for this figure
            if isempty(obj.hFig)
                obj.hFig = figure();
            end

            obj.clf();
            axes(obj.hFig);
            obj.hold('on');
        end

        function axis(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                axis(hAx, varargin{:});
            end
        end
        
        function val = axGet(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                if isempty(hAx)
                    val = [];
                else
                    val = get(hAx, varargin{:});
                end
            else
                val = [];
            end
        end
        
        function axSet(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                set(hAx, varargin{:});
            end
        end

        function clf(obj)
            if ~strcmp(obj.status, 'destroyed')
                clf(obj.hFig);
            end
        end

        function close(obj)
            %CLOSE Close the figure
            if ~strcmp(obj.status, 'destroyed')
                close(obj.hFig);
            end
        end

        function colorbar(obj, varargin)
            %COLORBAR Colorbar showing color scale
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                colorbar(hAx);
            end
        end
        
        function val = figGet(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                val = get(obj.hFig, varargin{:});
            else
                val = [];
            end
        end

        function figSet(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                set(obj.hFig, varargin{:});
            end
        end

        function grid(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                grid(hAx, varargin{:});
            end
        end

        function hold(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                hold(hAx, varargin{:});
            end
        end

        function plot(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                if isempty(hAx)
                    obj.toForeground();
                    plot(varargin{:});
                else
                    plot(hAx, varargin{:});
                end

                if obj.isMouseable
                    obj.setMouseable();
                end
            end
        end

        function setMouseable(obj, hFunClick)
            %SETMOUSEABLE Set the figure to be mouseable
            obj.isMouseable = true;
            if nargin == 2 && ~isempty(hFunClick)
                obj.hFunClick = hFunClick;
            end                

            if ~strcmp(obj.status, 'destroyed')
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
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                title(hAx, t, 'Interpreter', 'none', 'FontWeight', 'normal');
            end
        end

        function toForeground(obj)
            if ~strcmp(obj.status, 'destroyed')
                figure(obj.hFig);
            end
        end

        function hMenu = uimenu(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hMenu = uimenu(obj.hFig, varargin{:});
            else
                hMenu = [];
            end
        end

        function xlabel(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                xlabel(hAx, varargin{:});
            end
        end

        function ylabel(obj, varargin)
            if ~strcmp(obj.status, 'destroyed')
                hAx = get(obj.hFig, 'CurrentAxes');
                ylabel(hAx, varargin{:});
            end
        end
    end

    %% UTILITY METHODS
    methods (Hidden)
        function destroyFig(obj, hObject, ~)
            delete(hObject);
            obj.status = 'destroyed';
        end
    end

    %% GETTERS/SETTERS
    methods
        % figMetadata
        function fm = get.figMetadata(obj)
            if ~strcmp(obj.status, 'destroyed')
                fm = get(obj.hFig, 'UserData');
            else
                fm = [];
            end
        end
        function set.figMetadata(obj, fm)
            if ~strcmp(obj.status, 'destroyed')
                set(obj.hFig, 'UserData', fm);
            end
        end

        % figName
        function fn = get.figName(obj)
            if ~strcmp(obj.status, 'destroyed')
                fn = get(obj.hFig, 'Name');
            else
                fn = '';
            end
        end
        function set.figName(obj, figName)
            if ~strcmp(obj.status, 'destroyed')
                set(obj.hFig, 'Name', figName);
            end
        end

        % figPos
        function fp = get.figPos(obj)
            if ~strcmp(obj.status, 'destroyed')
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
            if ~strcmp(obj.status, 'destroyed')
                ft = get(obj.hFig, 'Tag');
            else
                ft = '';
            end
        end
        function set.figTag(obj, figTag)
            if ~strcmp(obj.status, 'destroyed')
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
            if ~strcmp(obj.status, 'destroyed')
                set(obj.hFig, 'KeyPressFcn', obj.hFunKey);
            end
        end

        % outerPosition
        function op = get.outerPosition(obj)
            if ~strcmp(obj.status, 'destroyed')
                op = get(obj.hFig, 'OuterPosition');
            else
                op = [];
            end
        end
        function set.outerPosition(obj, op)
            if ~strcmp(obj.status, 'destroyed')
                set(obj.hFig, 'OuterPosition', op);
            end
        end
    end
end

