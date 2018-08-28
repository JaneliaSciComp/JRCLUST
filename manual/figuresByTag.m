%--------------------------------------------------------------------------
function [hFig, figData] = figuresByTag(tag, hFig)
    % return figure handle(s) based on tags given
    % cache figure handles
    % [usage]
    % figuresByTag(tag)
    % figuresByTag(tag, hFig) %build cache

    if iscell(tag) % multiple tags
        hFig = nan(size(tag));
        for iFig = 1:numel(tag)
            try
                hFig(iFig) = findobj('Tag', tag{iFig});
            catch
            end
        end
    else
        figData = [];
        try
            hFig = findobj('Tag', tag, 'Type', 'figure');

            if isempty(hFig)
                hFig = createFigure(tag);

                try % set position if exists
                    S0 = get(0, 'UserData');
                    iFig = find(strcmp(S0.figTags, tag), 1, 'first');

                    if ~isempty(iFig)
                        set(hFig, 'OuterPosition', S0.figPositions{iFig});
                    end
                catch
                end
            else
                hFig = hFig(end); % get later one
            end

            if nargout > 1
                figData = get(hFig, 'UserData');
            end
        catch
            hFig = [];
            disperr_();
        end
    end
end % func
