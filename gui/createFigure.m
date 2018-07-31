%--------------------------------------------------------------------------
function hFig = createFigure(tag, figPos, figName, hasToolbar, hasMenubar)
    % or call external create_figure()
    if nargin < 2
        figPos = [];
    end

    if nargin < 3
        figName = '';
    end
    
    if nargin < 4
        hasToolbar = 0;
    end
    
    if nargin < 5
        hasMenubar = 0;
    end
    
    if isempty(tag)
        hFig = figure();
    else
        deleteMany(findobj('Tag', tag, 'Type', 'Figure'));
        hFig = figure('Tag', tag);
    end
    
    set(hFig, 'Name', figName, 'NumberTitle', 'off', 'Color', 'w');
    clf(hFig);

    set(hFig, 'UserData', []); % clear user data
    
    if ~hasToolbar
        set(hFig, 'ToolBar', 'none');
    else
        set(hFig, 'ToolBar', 'figure');
    end
    
    if ~hasMenubar
        set(hFig, 'MenuBar', 'none');
    else
        set(hFig, 'MenuBar', 'figure');
    end

    if ~isempty(figPos)
        resize_figure_(hFig, figPos);
    end
    clf(hFig);
end %func
