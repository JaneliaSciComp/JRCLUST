function hFigMap = doPlotFigMap(hFigMap, hClust, hCfg, selected)
    %DOPLOTFIGMAP Plot probe map
    iWaveforms = hClust.meanWfGlobal(:, :, selected(1));
    clusterVpp = squeeze(max(iWaveforms) - min(iWaveforms));
    vpp = repmat(clusterVpp(:)', [4, 1]);

    if isempty(hFigMap.figData)
        hFigMap.axes();
        [XData, YData] = getPadCoordinates(hCfg);
        hFigMap.addPatch('hPatch', XData, YData, vpp, 'EdgeColor', 'k', 'FaceColor', 'flat');
        hFigMap.axis([min(XData(:)), max(XData(:)), min(YData(:)), max(YData(:))]);
        % colormap jet;
        hFigMap.setMouseable();
        nSites = numel(hCfg.siteMap);
        hFigMap.addText('hText', hCfg.siteLoc(:, 1), hCfg.siteLoc(:, 2), ...
                        arrayfun(@(i) num2str(i), 1:nSites, 'UniformOutput', false), ...
                        'VerticalAlignment', 'bottom', ...
                        'HorizontalAlignment', 'left');
        hFigMap.xlabel('X Position (\mum)');
        hFigMap.ylabel('Y Position (\mum)');
    else
        hFigMap.plotSet('hPatch', 'CData', vpp);
    end

    hFigMap.title(sprintf('Max: %0.1f \\muVpp', max(clusterVpp)));
    hFigMap.caxis([0, max(clusterVpp)]);
end

%% LOCAL FUNCTIONS
function [XData, YData] = getPadCoordinates(hCfg)
    xOffsets = [0 0 1 1] * hCfg.probePad(2);
    yOffsets = [0 1 1 0] * hCfg.probePad(1);
    XData = bsxfun(@plus, hCfg.siteLoc(:, 1)', xOffsets(:));
    YData = bsxfun(@plus, hCfg.siteLoc(:, 2)', yOffsets(:));
end
