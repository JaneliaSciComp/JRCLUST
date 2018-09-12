%--------------------------------------------------------------------------
% 8/6/17 JJJ: Generalized to any figures previously querried.
function [hFig, S_fig] = get_fig_cache_(vcFig_tag)
    % clear before starting manual
    % Return from persistent cache
    % Create a new figure if Tag doesn't exist

    persistent S_fig_cache_
    % persistent hFigPos hFigMap hFigWav hFigTime hFigProj hFigWavCor hFigHist hFigIsi hFigCorr hFigRD hFig_traces hFig_preview
    if iscell(vcFig_tag)
        if nargout==1
            hFig = cellfun(@(vc)get_fig_cache_(vc), vcFig_tag, 'UniformOutput', 0);
        else
            [hFig, S_fig] = cellfun(@(vc)get_fig_cache_(vc), vcFig_tag, 'UniformOutput', 0);
        end
        return;
    end
    if isempty(S_fig_cache_)
        hFig = get_fig_(vcFig_tag);
        S_fig_cache_ = struct(vcFig_tag, hFig);
    else
        if isfield(S_fig_cache_, vcFig_tag)
            hFig = S_fig_cache_.(vcFig_tag);
            if isvalid_(hFig)
                hFig = S_fig_cache_.(vcFig_tag);
            else
                hFig = get_fig_(vcFig_tag);
                S_fig_cache_.(vcFig_tag) = hFig;
            end
        else
            hFig = get_fig_(vcFig_tag);
            S_fig_cache_.(vcFig_tag) = hFig;
        end
    end
    if nargout>1, S_fig = get(hFig, 'UserData'); end
end %func
