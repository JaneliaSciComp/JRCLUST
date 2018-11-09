%--------------------------------------------------------------------------
function hFig = create_figure_(figTag, figPos, figName, figToolbar, figMenubar)
    % create_figure_('FigPos', [0 0 .15 .5], ['Unit position; ', P.vcFile_prm], 1, 0);
    % or call external create_figure()
    if nargin < 2
        figPos = [];
    end
    if nargin < 3
        figName = '';
    end
    if nargin < 4
        figToolbar = false;
    end
    if nargin < 5
        figMenubar = false;
    end
%     if isempty(vcTag)
%         hFig = figure();
%     elseif ischar(vcTag)
%         hFig = figure_new_(vcTag);
%     else
%         hFig = vcTag;
%     end
%     set(hFig, 'Name', vcName, 'NumberTitle', 'off', 'Color', 'w');
%     clf(hFig);
%     set(hFig, 'UserData', []); %empty out the user data
%     if ~fToolbar
%         set(hFig, 'ToolBar', 'none');
%     else
%         set(hFig, 'ToolBar', 'figure');
%     end
%     if ~fMenubar
%         set(hFig, 'MenuBar', 'none');
%     else
%         set(hFig, 'MenuBar', 'figure');
%     end
% 
%     if ~isempty(vrPos), resize_figure_(hFig, vrPos); end
%     clf(hFig);

    hFig = jrclust.views.Figure(figTag, figPos, figName, figToolbar, figMenubar);
end %func
