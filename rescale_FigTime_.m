%--------------------------------------------------------------------------
function rescale_FigTime_(event, S0, P)
    % rescale_FigTime_(event, S0, P)
    % rescale_FigTime_(maxAmp, S0, P)
    if nargin<2, S0 = []; end
    if nargin<3, P = []; end

    if isempty(S0), S0 = get0_(); end
    if isempty(P), P = S0.P; end
    S_clu = S0.S_clu;

    [S_fig, maxAmp_prev] = set_fig_maxAmp_('FigTime', event);
    ylim_(S_fig.hAx, [0, 1] * S_fig.maxAmp);
    imrect_set_(S_fig.hRect, [], [0, S_fig.maxAmp]);
    iSite = S_clu.clusterSites(S0.iCluCopy);

    % switch lower(P.displayFeature)
    %     case {'vpp', 'vmin'} %voltage feature
    %         vcYlabel = sprintf('Site %d (\\mu%s)', iSite, P.displayFeature);
    %     otherwise %other feature options
    %         vcYlabel = sprintf('Site %d (%s)', iSite, P.displayFeature);
    % end
    % ylabel(S_fig.hAx, vcYlabel);

end %func
