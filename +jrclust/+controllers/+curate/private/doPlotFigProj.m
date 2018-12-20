function hFigProj = doPlotFigProj(hFigProj, hClust, hCfg, selected, maxAmp)
    %DOPLOTFIGPROJ Plot feature projection figure
    hFigProj.rmPlot('hSelect'); % clear select polygon
    hFigProj.hidePlot('background'); % clear select polygon

    iSite = hClust.clusterSites(selected(1));
    
    % limit the number of sites to display in the feature projection view
    nSites = min(hCfg.nSitesFigProj, size(hCfg.siteNeighbors, 1)); % by request

    % center sites around cluster center site
    if nSites < size(hCfg.siteNeighbors, 1)
        sitesToShow = iSite:iSite + nSites - 1;
        if sitesToShow(end) > max(hCfg.siteMap) % correct for overshooting
            sitesToShow = sitesToShow - max(sitesToShow) + max(hCfg.siteMap);
        end
    else
        sitesToShow = sort(hCfg.siteNeighbors(:, iSite), 'ascend');
    end

    if strcmp(hCfg.dispFeature, 'vpp')
        xLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        yLabel = 'Site # (%0.0f \\muV_{min})';    
    elseif ismember(hCfg.dispFeature, {'kilosort', 'pca', 'gpca', 'ppca'})
        xLabel = sprintf('Site # (PC %d)', 1); % fix this in post
        yLabel = sprintf('Site # (PC %d)', 2);
    else
        xLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', hCfg.dispFeature, hCfg.dispFeature, hCfg.dispFeature);
        yLabel = sprintf('Site # (%%0.0f %s)', hCfg.dispFeature);
    end
    figTitle = '[H]elp; [S]plit; [B]ackground; (Sft)[Up/Down]:Scale; [Left/Right]:Sites; [M]erge; [F]eature';

    if isempty(hFigProj.figData)
        hFigProj.axes();
        hFigProj.axSet('Position', [.1 .1 .85 .85], 'XLimMode', 'manual', 'YLimMode', 'manual');

        hFigProj.addLine('background',  nan, nan, 'Color', hCfg.mrColor_proj(1, :));
        hFigProj.addLine('foreground',  nan, nan, 'Color', hCfg.mrColor_proj(2, :)); % placeholder
        hFigProj.addLine('foreground2', nan, nan, 'Color', hCfg.mrColor_proj(3, :)); % placeholder
        
        plotStyle = {'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'};
        hFigProj.plotSet('background', plotStyle{:});
        hFigProj.plotSet('foreground', plotStyle{:});
        hFigProj.plotSet('foreground2', plotStyle{:});

        % plot boundary
        hFigProj.addTable('hTable', [0, nSites], '-', 'Color', [.5 .5 .5]);
        hFigProj.addDiag('hDiag', [0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5);
        hFigProj.setHideOnDrag('background');
        hFigProj.figData.isPlotted = true;
    end

    dispFeatures = getDispFeatures(hClust, hCfg, sitesToShow, selected);
    bgYData = dispFeatures.bgY;
    bgXData = dispFeatures.bgX;
    fgYData = dispFeatures.fgY;
    fgXData = dispFeatures.fgX;
    fg2YData = dispFeatures.fg2Y;
    fg2XData = dispFeatures.fg2X;

%     if ~isfield(hFigProj.figData, 'viSites_show')
%         hFigProj.figData.viSites_show = [];
%     end

    %if ~equal_vr_(hFigProj.figData.viSites_show, hCfg.viSites_show) || ~equal_vr_(hFigProj.figData.dispFeature, hCfg.viSites_show)
    % plot background spikes
    plotFeatures(hFigProj, 'background', bgYData, bgXData, maxAmp, hCfg);
    %end

    % plot foreground spikes
    plotFeatures(hFigProj, 'foreground', fgYData, fgXData, maxAmp, hCfg);
    % plot secondary foreground spikes
    if numel(selected) == 2
        plotFeatures(hFigProj, 'foreground2', fg2YData, fg2XData, maxAmp, hCfg);
        figTitle = sprintf('Clu%d (black), Clu%d (red); %s', selected(1), selected(2), figTitle);
    else % or hide the plot
        hFigProj.hidePlot('background');
        figTitle = sprintf('Clu%d (black); %s', selected(1), figTitle);
    end

    % Annotate axes
    hFigProj.axis([0 nSites 0 nSites]);
    hFigProj.axSet('XTick', 0.5:1:nSites, 'YTick', 0.5:1:nSites, ...
                   'XTickLabel', sitesToShow, 'YTickLabel', sitesToShow, ...
                   'Box', 'off');
    hFigProj.xlabel(sprintf(xLabel, maxAmp));
    hFigProj.ylabel(sprintf(yLabel, maxAmp));
    hFigProj.title(figTitle);

    hFigProj.figData.csHelp = {'[D]raw polygon', ...
                    '[S]plit cluster', ...
                    '(shift)+Up/Down: change scale', ...
                    '[R]eset scale', ...
                    'Zoom: mouse wheel', ...
                    'Drag while pressing wheel: pan'};
end

%% LOCAL FUNCTIONS
function plotFeatures(hFigProj, plotKey, featY, featX, boundScale, hCfg)
%     switch lower(P.dispFeature)
%         case 'vpp'
%             bounds = maxAmp*[0 1];
% 
%         otherwise
%             % round up to nearest 50 on either side of 0
%             bounds = maxAmp*[-1 1];
%     end

    bounds = boundScale*[0 1];
    [XData, YData] = ampToProj(featY, featX, bounds, hCfg.nSiteDir, hCfg);
    hFigProj.updatePlot(plotKey, XData, YData);
end
