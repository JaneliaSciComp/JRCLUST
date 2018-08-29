%--------------------------------------------------------------------------
% 122917 JJJ: modified
function figData = rescaleFigProj(event, hFig, figData, S0)

    if nargin < 2
        hFig = [];
    end
    if nargin < 3
        figData = [];
    end
    if nargin < 4 || isempty(S0)
        S0 = get(0, 'UserData');
    end

    if isempty(hFig) || isempty(figData)
        [hFig, figData] = getCachedFig('FigProj');
    end

    P = S0.P;

    if isnumeric(event)
        figData.maxAmp = event;
    else
        figData.maxAmp = change_amp_(event, figData.maxAmp);
    end

    maxAmp = figData.maxAmp;

    hPlots = [figData.hPlotBG, figData.hPlotFG, figData.hPlotFG2];
    if isempty(S0.secondarySelectedCluster)
        hPlots(end) = [];
    end

    % rescaleProj_(hPlots, figData.maxAmp, S0.P);
    clearPlots(figData.hPlotFG2); % replacement for update_plot2_proj_
    for iPlot = 1:numel(hPlots)
        hPlot = hPlots(iPlot);
        plotData = get(hPlot, 'UserData');

        if isfield(plotData, 'hPoly')
            delete(plotData.hPoly);
        end

        if isempty(plotData)
            continue;
        end

        switch lower(P.displayFeature)
            case 'kilosort'
                % round to nearest 50 on either side of 0
                bounds = round(max(abs([plotData.mrMin(:) ; plotData.mrMax(:)]))/50, 0)*[-50, 50];
                maxPair = [];

            otherwise % vpp et al.
                bounds = [0 maxAmp];
                maxPair = P.maxSite_show;
        end

        [vrX, vrY, viPlot, ~] = featuresToSiteGrid(plotData.mrMax, plotData.mrMin, bounds, maxPair);

        plotData = struct_add_(plotData, viPlot, vrX, vrY, maxAmp);
        updatePlot(hPlot, vrX, vrY, plotData);
    end

    switch lower(P.displayFeature)
        case {'vpp', 'vmin', 'vmax'}
            figData.vcXLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
            figData.vcYLabel = 'Site # (%0.0f \\muV_{min})';

        case 'kilosort'
            figData.vcXLabel = sprintf('Site # (%%0.0f KS PC %d)', S0.pcPair(1));
            figData.vcYLabel = sprintf('Site # (%%0.0f KS PC %d)', S0.pcPair(2));

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
