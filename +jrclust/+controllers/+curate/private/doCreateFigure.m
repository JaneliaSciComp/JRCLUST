function hFig = doCreateFigure(figTag, figPos, figName, figToolbar, figMenubar)
    %DOCREATEFIGURE Create a Figure object
    if nargin < 2
        figPos = [];
    end
    if nargin < 3
        figName = '';
    end
    if nargin < 4
        figToolbar = 0;
    end
    if nargin < 5
        figMenubar = 0;
    end

    hFig = jrclust.views.Figure(figTag, figPos, figName, figToolbar, figMenubar);
end
