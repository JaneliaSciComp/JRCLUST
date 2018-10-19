%--------------------------------------------------------------------------
function [S_fig, maxAmp_prev, hFig] = set_fig_maxAmp_(vcFig, event)
    [hFig, S_fig] = get_fig_cache_(vcFig);
    if isempty(S_fig)
        P = get0_('P');
        S_fig.maxAmp = P.maxAmp;
    end
    maxAmp_prev = S_fig.maxAmp;
    if isnumeric(event)
        S_fig.maxAmp = event;
    else
        S_fig.maxAmp = change_amp_(event, maxAmp_prev);
    end
    set(hFig, 'UserData', S_fig);
end
