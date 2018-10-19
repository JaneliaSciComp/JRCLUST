function MainFig = mouse_figure(MainFig, axs, hFunClick)
%MOUSE_FIGURE           mouse-friendly figure
%
% MOUSE_FIGURE() creates a figure (or modifies an existing one) that allows
% zooming with the scroll wheel and panning with mouse clicks, *without*
% first selecting the ZOOM or PAN tools from the toolbar. Moreover, zooming
% occurs to and from the point the mouse currently hovers over, instead of
% to and from the less intuitive "CameraPosition". 
%
%         Scroll: zoom in/out
%     Left click: pan
%   Double click: reset view to default view
%    Right click: set new default view
%
% LIMITATIONS: This function (re-)efines several functions in the figure 
% (WindowScrollWheelFcn, WindowButtonDownFcn, WindowButtonUpFcn and 
% WindowButtonMotionFcn), so if you have any of these functions already
% defined they will get overwritten. Also, MOUSE_FIGURE() only works 
% properly for 2-D plots. As such, it should only be used for simple, 
% first-order plots intended for "quick-n-dirty" data exploration.
%
% EXAMPLE:
%
%   mouse_figure;
%   x = linspace(-1, 1, 10000);
%   y = sin(1./x);
%   plot(x, y) 
%   
% See also figure, axes, zoom, pan.


% Please report bugs and inquiries to: 
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com    (personal)
%              oldenhuis@luxspace.lu  (professional)
% Affiliation: LuxSpace sàrl
% Licence    : BSD


% If you find this work useful, please consider a donation:
% https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6G3S5UYM7HJ3N

    
    % initialize
    status = '';  previous_point = [];
    
    % initialize axes
    if (nargin == 0) || ~ishandle(MainFig)
        MainFig = figure;  axs = gca;
    elseif nargin==1
        axs = get(MainFig, 'currentaxes');
    end
    if nargin<3, hFunClick = []; end
    
    % only works properly for 2D plots
    if ~is2D(axs) % is2D might disappear in a future release...
        error('mouse_figure:plot3D_not_supported', ...
              'MOUSE_FIGURE() only works for 2-D plots.');
    end
        
    % get original limits
    original_xlim = get(axs, 'xlim');
    original_ylim = get(axs, 'ylim');
    
    % define zooming with scrollwheel, and panning with mouseclicks
    set(MainFig, ...
        'WindowScrollWheelFcn' , @scroll_zoom,...
        'WindowButtonDownFcn'  , @pan_click,...
        'WindowButtonUpFcn'    , @pan_release,...
        'WindowButtonMotionFcn', @pan_motion);    
%     set(axs, 'buttonDownFcn', @button_help); %JJJ @TODO
%     
%     % JJJ
%     function button_help(hObject, event)
%         switch lower(event.Key)
%             case 'h'
%             csHelp = { ...
%                 'h: show help', ...
%                 '------------------', ...
%                 'middle-click: drag', ...
%                 'Wheel: zoom in/out'}; 
%             uiwait(msgbox(csHelp, 'modal'));
%         end %switch
%     end %func

    % zoom in to the current point with the mouse wheel
    function scroll_zoom(varargin)
%         hMsg = msgbox_open('Zooming...');
        % double check if these axes are indeed the current axes
        if get(MainFig, 'currentaxes') ~= axs, return, end
        % get the amount of scolls
        scrolls = varargin{2}.VerticalScrollCount;
        % get the axes' x- and y-limits
        xlim = get(axs, 'xlim');  ylim = get(axs, 'ylim');
        % get the current camera position, and save the [z]-value
        cam_pos_Z = get(axs, 'cameraposition');  cam_pos_Z = cam_pos_Z(3);
        % get the current point
        old_position = get(axs, 'CurrentPoint'); old_position(1,3) = cam_pos_Z;
        % calculate zoom factor
        if scrolls < 0
            zoomfactor = sqrt(2); %1 + 1/4; %zoom speed
        else
            zoomfactor = 1/sqrt(2); %1 - 1/4; %zoom speed
        end
        % adjust camera position
%         set(axs, 'Visible', 'off');
        set(axs, 'cameratarget', [old_position(1, 1:2), 0],...
            'cameraposition', old_position(1, 1:3));
        % adjust the camera view angle (equal to zooming in)
        camzoom(zoomfactor);
        keyFig = get(MainFig, 'CurrentCharacter');

        [zoomfactor_x, zoomfactor_y] = deal(zoomfactor);
        if lower(keyFig) == 'x', zoomfactor_y = 1; end
        if lower(keyFig) == 'y', zoomfactor_x = 1; end
        
        % zooming with the camera has the side-effect of
        % NOT adjusting the axes limits. We have to correct for this:
        x_lim1 = (old_position(1,1) - min(xlim))/zoomfactor_x;
        x_lim2 = (max(xlim) - old_position(1,1))/zoomfactor_x;
        xlim   = [old_position(1,1) - x_lim1, old_position(1,1) + x_lim2];
        y_lim1 = (old_position(1,2) - min(ylim))/zoomfactor_y;
        y_lim2 = (max(ylim) - old_position(1,2))/zoomfactor_y;
        ylim   = [old_position(1,2) - y_lim1, old_position(1,2) + y_lim2];
        set(axs, 'xlim', xlim), set(axs, 'ylim', ylim)
        % set new camera position
        new_position = get(axs, 'CurrentPoint');
        old_camera_target =  get(axs, 'CameraTarget');
        old_camera_target(3) = cam_pos_Z;
        new_camera_position = old_camera_target - ...
            (new_position(1,1:3) - old_camera_target(1,1:3));
        % adjust camera target and position
        set(axs, 'cameraposition', new_camera_position(1, 1:3),...
            'cameratarget', [new_camera_position(1, 1:2), 0]);
        % we also have to re-set the axes to stretch-to-fill mode
        set(axs, 'cameraviewanglemode', 'auto',...
            'camerapositionmode', 'auto',...
            'cameratargetmode', 'auto');
%         msgbox_close(hMsg);
%         set(axs, 'Visible', 'on');
    end % scroll_zoom
    
    % pan upon mouse click
    function pan_click(varargin)
        % double check if these axes are indeed the current axes
        if get(MainFig, 'currentaxes') ~= axs, return, end
        % perform appropriate action
        vcButton = lower(get(MainFig, 'selectiontype')); %JJJ
        switch vcButton
            % start panning on left click
            case {'normal' , 'alt'}
                xyPoint = get(axs, 'CurrentPoint');
                if ~isempty(hFunClick)
                    hFunClick(xyPoint([1,3]), vcButton); 
                end
            case 'extend'
                status = 'down';
                % hide objects
                hide_drag_(MainFig);                
                previous_point = get(axs, 'CurrentPoint');                              
        end
    end
    
    % release mouse button
    function pan_release(varargin)
        show_drag_(MainFig);
%         set(axs, 'Visible', 'on');
        % double check if these axes are indeed the current axes
        if get(MainFig, 'currentaxes') ~= axs, return, end
        % just reset status
        status = '';
    end
    
    % move the mouse (with button clicked)
    function pan_motion(varargin)
        % double check if these axes are indeed the current axes
        if get(MainFig, 'currentaxes') ~= axs, return, end
        % return if there isn't a previous point
        if isempty(previous_point), return, end  
        % return if mouse hasn't been clicked
        if isempty(status), return, end  
        % get current location (in pixels)
        current_point = get(axs, 'CurrentPoint');
        % get current XY-limits
        xlim = get(axs, 'xlim');  ylim = get(axs, 'ylim');     
        % find change in position
        delta_points = current_point - previous_point;  
        % adjust limits
        new_xlim = xlim - delta_points(1); 
        new_ylim = ylim - delta_points(3); 
        % set new limits
        set(axs, 'Xlim', new_xlim); set(axs, 'Ylim', new_ylim);           
        % save new position
        previous_point = get(axs, 'CurrentPoint');
    end 
    
end


%--------------------------------------------------------------------------
function hide_drag_(hFig)
global mouse_figure_hidden %remember hidden obj
% S0 = get(0, 'UserData');
S_fig = get(hFig, 'UserData');
if ~isfield(S_fig, 'cvhHide_mouse'), return; end
try 
%     iFig = find(hFig == S0.vhFig_mouse);
%     if isempty(iFig), return; end
    vhHide = S_fig.cvhHide_mouse;
    if strcmpi(get(vhHide{1}(1), 'Visible'), 'off'), return; end
    
    toggleVisible_(vhHide, 0); %hide
    mouse_figure_hidden = vhHide; %hidden object
catch
    ;
end
end %func


%--------------------------------------------------------------------------
function show_drag_(hFig)
global mouse_figure_hidden
if isempty(mouse_figure_hidden), return; end 
try 
    toggleVisible_(mouse_figure_hidden, 1); %hide
    mouse_figure_hidden = [];
catch
    ;
end
end %func


%--------------------------------------------------------------------------
function vlVisible = toggleVisible_(vhPlot, fVisible)
if isempty(vhPlot), return; end

if iscell(vhPlot)
    cvhPlot = vhPlot;
    if nargin<2
        cellfun(@(vhPlot)toggleVisible_(vhPlot), cvhPlot);
    else
        cellfun(@(vhPlot)toggleVisible_(vhPlot, fVisible), cvhPlot);
    end
    return;
end
try
    if nargin==1
        vlVisible = false(size(vhPlot));
        % toggle visibility
        for iH=1:numel(vhPlot)
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
end %func