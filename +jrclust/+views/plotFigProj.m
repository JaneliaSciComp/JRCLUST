function hFigProj = plotFigProj(hFigProj, hClust, sitesToShow, selected, boundScale)
    %PLOTFIGPROJ Plot feature projection figure
    hCfg = hClust.hCfg;

    hFigProj.clearPlot('foreground2'); % clear secondary cluster spikes

    if strcmp(hCfg.dispFeature, 'vpp')
        XLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        YLabel = 'Site # (%0.0f \\muV_{min})';
    elseif ismember(hCfg.dispFeature, {'template', 'pca', 'ppca'})
        XLabel = sprintf('Site # (PC %d)', hCfg.pcPair(1));
        YLabel = sprintf('Site # (PC %d)', hCfg.pcPair(2));
    else
        XLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', hCfg.dispFeature, hCfg.dispFeature, hCfg.dispFeature);
        YLabel = sprintf('Site # (%%0.0f %s)', hCfg.dispFeature);
    end
    figTitle = '[H]elp; [S]plit; [B]ackground; (Sft)[Up/Down]:Scale; [Left/Right]:Sites; [M]erge; [F]eature';

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

        % save for later
        hFigProj.figData.boundScale = boundScale;
    end

    dispFeatures = getFigProjFeatures(hClust, sitesToShow, selected);
    bgYData = dispFeatures.bgYData;
    bgXData = dispFeatures.bgXData;
    fgYData = dispFeatures.fgYData;
    fgXData = dispFeatures.fgXData;
    fg2YData = dispFeatures.fg2YData;
    fg2XData = dispFeatures.fg2XData;

    % save these for autoscaling
    hFigProj.figData.dispFeatures = dispFeatures;

    % plot background spikes
    plotFeatures(hFigProj, 'background', bgYData, bgXData, boundScale, hCfg);

    % plot foreground spikes
    assigns = plotFeatures(hFigProj, 'foreground', fgYData, fgXData, boundScale, hCfg);
    hFigProj.figData.dispFeatures.fgSpikes = hFigProj.figData.dispFeatures.fgSpikes(assigns);

    % plot secondary foreground spikes
    if numel(selected) == 2
        plotFeatures(hFigProj, 'foreground2', fg2YData, fg2XData, boundScale, hCfg);
        figTitle = sprintf('Clu%d (black), Clu%d (red); %s', selected(1), selected(2), figTitle);
    else % or hide the plot
        hFigProj.clearPlot('foreground2');
        figTitle = sprintf('Clu%d (black); %s', selected(1), figTitle);
    end

    % Annotate axes
    hFigProj.axApply('default', @axis, [0 nSites 0 nSites]);
    hFigProj.axApply('default', @set, 'XTick', 0.5:1:nSites, 'YTick', 0.5:1:nSites, ...
                     'XTickLabel', sitesToShow, 'YTickLabel', sitesToShow, ...
                    'Box', 'off');
    hFigProj.axApply('default', @xlabel, sprintf(XLabel, boundScale));
    hFigProj.axApply('default', @ylabel, sprintf(YLabel, boundScale));
    hFigProj.axApply('default', @title, figTitle);

    hFigProj.figData.helpText = {'[D]raw polygon', ...
                    '[S]plit cluster', ...
                    '(shift)+Up/Down: change scale', ...
                    '[R]eset scale', ...
                    'Zoom: mouse wheel', ...
                    'Drag while pressing wheel: pan'};
end

%% LOCAL FUNCTIONS
function assigns = plotFeatures(hFigProj, plotKey, featY, featX, boundScale, hCfg)
    %PLOTFEATURES Plot features in a grid
    if strcmp(hCfg.dispFeature, 'vpp')
        bounds = boundScale*[0 1];
    else
        bounds = boundScale*[-1 1];
    end

    [XData, YData, assigns] = ampToProj(featY, featX, bounds, hCfg.nSiteDir, hCfg);
    hFigProj.updatePlot(plotKey, XData, YData);
end