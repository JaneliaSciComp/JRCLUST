%--------------------------------------------------------------------------
% 122917 JJJ: modified
function S_fig = rescale_FigProj_(event, hFig, S_fig, S0)
    % S_fig = rescale_FigProj_(event, hFig, S_fig, S0)
    % S_fig = rescale_FigProj_(maxAmp)

    if nargin<2, hFig = []; end
    if nargin<3, S_fig = []; end
    if nargin<4, S0 = []; end
    if isempty(hFig) || isempty(S_fig), [hFig, S_fig] = get_fig_cache_('FigProj'); end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P;

    if isnumeric(event)
        S_fig.maxAmp = event;
    else
        S_fig.maxAmp = change_amp_(event, S_fig.maxAmp);
    end
    vhPlot = [S_fig.hPlot0, S_fig.hPlot1, S_fig.hPlot2];
    if isempty(S0.iCluPaste), vhPlot(end) = []; end
    rescaleProj_(vhPlot, S_fig.maxAmp, S0.P);
    switch lower(P.vcFet_show)
        case {'vpp', 'vmin', 'vmax'}
        S_fig.vcXLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        S_fig.vcYLabel = 'Site # (%0.0f \\muV_{min})';
        otherwise
        S_fig.vcXLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', P.vcFet_show, P.vcFet_show, P.vcFet_show);
        S_fig.vcYLabel = sprintf('Site # (%%0.0f %s)', P.vcFet_show);
    end
    xlabel(S_fig.hAx, sprintf(S_fig.vcXLabel, S_fig.maxAmp));
    ylabel(S_fig.hAx, sprintf(S_fig.vcYLabel, S_fig.maxAmp));
    if nargout==0
        set(hFig, 'UserData', S_fig);
    end
end
