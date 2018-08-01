%--------------------------------------------------------------------------
% 8/6/17 JJJ: Generalized to any figures previously queried.
function [hFig, figData] = getCachedFig(tag)
    % clear before starting manual
    % Return from persistent cache
    % Create a new figure if Tag doesn't exist

    persistent figCache;

    if iscell(tag)
        if nargout == 1
            hFig = cellfun(@(c) getCachedFig(c), tag, 'UniformOutput', 0);
        else
            [hFig, figData] = cellfun(@(c) getCachedFig(c), tag, 'UniformOutput', 0);
        end
    else
        if isempty(figCache) % initialize figure cache
            hFig = figureByTag(tag);
            figCache = struct(tag, hFig);
        else
            if isfield(figCache, tag)
                hFig = figCache.(tag);

                if tryIsValid(hFig) % check cached figure handle still valid
                    hFig = figCache.(tag);
                else % not valid, recache
                    hFig = figureByTag(tag);
                    figCache.(tag) = hFig;
                end
            else
                hFig = figureByTag(tag);
                figCache.(tag) = hFig;
            end
        end

        if nargout > 1
            figData = get(hFig, 'UserData');
        end
    end
end %func
