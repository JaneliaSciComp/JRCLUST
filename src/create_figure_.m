%--------------------------------------------------------------------------
function hFig = create_figure_(vcTag, vrPos, vcName, fToolbar, fMenubar)
    % or call external create_figure()
    if nargin<2, vrPos = []; end
    if nargin<3, vcName = ''; end
    if nargin<4, fToolbar = 0; end
    if nargin<5, fMenubar = 0; end
    if isempty(vcTag)
        hFig = figure();
    elseif ischar(vcTag)
        hFig = figure_new_(vcTag);
    else
        hFig = vcTag;
    end
    set(hFig, 'Name', vcName, 'NumberTitle', 'off', 'Color', 'w');
    clf(hFig);
    set(hFig, 'UserData', []); %empty out the user data
    if ~fToolbar
        set(hFig, 'ToolBar', 'none');
    else
        set(hFig, 'ToolBar', 'figure');
    end
    if ~fMenubar
        set(hFig, 'MenuBar', 'none');
    else
        set(hFig, 'MenuBar', 'figure');
    end

    if ~isempty(vrPos), resize_figure_(hFig, vrPos); end
    clf(hFig);
end %func
