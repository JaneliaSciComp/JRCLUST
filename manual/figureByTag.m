%--------------------------------------------------------------------------
function [hFig, S_fig] = figureByTag(tag, hFig)
    % return figure handle based on the tag
    % cache figure handles
    % [usage]
    % figureByTag(tag)
    % figureByTag(tag, hFig) %build cache

    % multiple tags requested
    if iscell(tag), hFig = get_fig_all_(tag); return; end

    S_fig = [];
    try
        hFig = findobj('Tag', tag, 'Type', 'figure');
        if isempty(hFig)
            hFig = createFigure(tag);
            try %set position if exists
                S0 = get(0, 'UserData');
                iFig = find(strcmp(S0.figTags, tag), 1, 'first');
                if ~isempty(iFig)
                    set(hFig, 'OuterPosition', S0.figPositions{iFig});
                end
            catch
                %             disperr_();
            end
        else
            hFig=hFig(end); %get later one
        end
        if nargout>1, S_fig = get(hFig, 'UserData'); end
    catch
        hFig = [];
        disperr_();
    end
end %end
