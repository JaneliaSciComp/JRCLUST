function hFigProj = plotFigProj(hFigProj, hClust, sitesToShow, selected, boundScale, channel_idx, doAutoscale)
    %PLOTFIGPROJ Plot feature projection figure
    if nargin < 7
        doAutoscale = 1;
    end

    hCfg = hClust.hCfg;

    hFigProj.clearPlot('foreground2'); % clear secondary cluster spikes

    if strcmp(hCfg.dispFeature, 'vpp')
        XLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        YLabel = 'Site # (%0.0f \\muV_{min})';
    elseif ismember(hCfg.dispFeature, {'template', 'pca', 'ppca'})
        XLabel = sprintf('Site # (PC %d) (a.u.)', hCfg.pcPair(1));
        YLabel = sprintf('Site # (PC %d) (a.u.)', hCfg.pcPair(2));
    else
        XLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', hCfg.dispFeature, hCfg.dispFeature, hCfg.dispFeature);
        YLabel = sprintf('Site # (%%0.0f %s)', hCfg.dispFeature);
    end

    nSites = numel(sitesToShow);
    if ~hFigProj.hasAxes('default')
        hFigProj.addAxes('default');
        hFigProj.axApply('default', @set, 'Position', [.1 .1 .85 .85], 'XLimMode', 'manual', 'YLimMode', 'manual');

        hFigProj.addPlot('background', @line, nan, nan, 'Color', hCfg.colorMap(1, :));
        hFigProj.addPlot('foreground', @line, nan, nan, 'Color', hCfg.colorMap(2, :)); % placeholder
        hFigProj.addPlot('foreground2', @line,  nan, nan, 'Color', hCfg.colorMap(3, :)); % placeholder

        plotStyle = {'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'};
        hFigProj.plotApply('background', @set, plotStyle{:});
        hFigProj.plotApply('foreground', @set, plotStyle{:});
        hFigProj.plotApply('foreground2', @set, plotStyle{:});

        % plot boundary
        hFigProj.addTable('hTable', [0, nSites], '-', 'Color', [.5 .5 .5]);
        hFigProj.addDiag('hDiag', [0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5);
        hFigProj.setHideOnDrag('background');
    end

    dispFeatures = getFigProjFeatures(hClust, sitesToShow, selected);
    bgYData = dispFeatures.bgYData;
    bgXData = dispFeatures.bgXData;
    fgYData = dispFeatures.fgYData;
    fgXData = dispFeatures.fgXData;
    fg2YData = dispFeatures.fg2YData;
    fg2XData = dispFeatures.fg2XData;

    if doAutoscale
        autoscalePct = hClust.hCfg.getOr('autoscalePct', 99.5)/100;

        if numel(selected) == 1
            projData = {fgYData, fgXData};
        else
            projData = {fgYData, fgXData, fg2YData, fg2XData};
        end

        if ~all(cellfun(@isempty, projData))
            projData = cellfun(@(x) x(~isnan(x)), projData, 'UniformOutput', 0);
            boundScale = max(cellfun(@(x) quantile(abs(x(:)), autoscalePct), projData));
        end
    end

    % save scales for later
    hFigProj.figData.initialScale = boundScale;
    hFigProj.figData.boundScale = boundScale;

    % plot background spikes
    plotFeatures(hFigProj, 'background', bgYData, bgXData, boundScale, hCfg);

    % plot foreground spikes
    plotFeatures(hFigProj, 'foreground', fgYData, fgXData, boundScale, hCfg);

    % plot secondary foreground spikes
    if numel(selected) == 2
        plotFeatures(hFigProj, 'foreground2', fg2YData, fg2XData, boundScale, hCfg);
        figTitle = sprintf('Unit %d (black), Unit %d (red); (press [H] for help)', selected(1), selected(2));
    else % or hide the plot
        hFigProj.clearPlot('foreground2');
        figTitle = sprintf('Unit %d (black); (press [H] for help)', selected(1));
        hFigProj.figData.foreground2.XData = nan;
        hFigProj.figData.foreground2.YData = nan;
    end

    % Annotate axes
    hFigProj.axApply('default', @axis, [0 nSites 0 nSites]);
    hFigProj.axApply('default', @set, 'XTick', 0.5:1:nSites, 'YTick', 0.5:1:nSites, ...
                     'XTickLabel', channel_idx(sitesToShow), 'YTickLabel', channel_idx(sitesToShow), ...
                    'Box', 'off');
    hFigProj.axApply('default', @xlabel, sprintf(XLabel, boundScale));
    hFigProj.axApply('default', @ylabel, sprintf(YLabel, boundScale));
    hFigProj.axApply('default', @title, figTitle);
end

%% LOCAL FUNCTIONS
function plotFeatures(hFigProj, plotKey, featY, featX, boundScale, hCfg)
    %PLOTFEATURES Plot features in a grid
    if strcmp(hCfg.dispFeature, 'vpp')
        bounds = boundScale*[0 1];
    else
        bounds = boundScale*[-1 1];
    end

    [XData, YData] = ampToProj(featY, featX, bounds, hCfg.nSiteDir, hCfg);

    hFigProj.figData.(plotKey).XData = XData;
    hFigProj.figData.(plotKey).YData = YData;

    % subset features if plotting foreground
    if strcmp(plotKey, 'foreground')
        nSpikes = size(YData, 1);
        subset = jrclust.utils.subsample(1:nSpikes, hCfg.nSpikesFigProj);
        XData = XData(subset, :);
        YData = YData(subset, :);
        hFigProj.figData.foreground.subset = subset;
    end
    hFigProj.updatePlot(plotKey, XData(:), YData(:));
end
