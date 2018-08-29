%--------------------------------------------------------------------------
% 122917 JJJ: modified
function figData = rescaleFigProj(event, hFig, figData, S0)
    % figData = rescaleFigProj(event, hFig, figData, S0)
    % figData = rescaleFigProj(maxAmp)

    if nargin < 2
        hFig = [];
    end
    if nargin < 3
        figData = [];
    end
    if nargin<4, S0 = []; end
    if isempty(hFig) || isempty(figData), [hFig, figData] = getCachedFig('FigProj'); end
    if isempty(S0), S0 = get(0, 'UserData'); end
    P = S0.P;

    if isnumeric(event)
        figData.maxAmp = event;
    else
        figData.maxAmp = change_amp_(event, figData.maxAmp);
    end % hPlotBG: background; hPlotFG: foreground

    vhPlot = [figData.hPlotBG, figData.hPlotFG, figData.hPlotFG2];

    if isempty(S0.secondarySelectedCluster)
        vhPlot(end) = [];
    end

    rescaleProj_(vhPlot, figData.maxAmp, S0.P);

    switch lower(P.displayFeature)
        case {'vpp', 'vmin', 'vmax'}
            figData.vcXLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
            figData.vcYLabel = 'Site # (%0.0f \\muV_{min})';

        case 'kilosort'
            figData.vcXLabel = sprintf('Site # (%%0.0f KS PC %d)', S0.kspc(1));
            figData.vcYLabel = sprintf('Site # (%%0.0f KS PC %d)', S0.kspc(2));

        otherwise
            figData.vcXLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', P.displayFeature, P.displayFeature, P.displayFeature);
            figData.vcYLabel = sprintf('Site # (%%0.0f %s)', P.displayFeature);
    end

    xlabel(figData.hAx, sprintf(figData.vcXLabel, figData.maxAmp));
    ylabel(figData.hAx, sprintf(figData.vcYLabel, figData.maxAmp));

    if nargout == 0
        set(hFig, 'UserData', figData);
    end
end
