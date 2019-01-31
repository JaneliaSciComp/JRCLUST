function hFigMap = doPlotFigMap(hFigMap, hClust, hCfg, selected)
    %DOPLOTFIGMAP Plot probe map
    iWaveforms = hClust.meanWfGlobal(:, :, selected(1));
    clusterVpp = squeeze(max(iWaveforms) - min(iWaveforms));
    vpp = repmat(clusterVpp(:)', [4, 1]);

    if ~hFigMap.hasAxes('default')
        hFigMap.addAxes('default');
        [XData, YData] = getPadCoordinates(hCfg);
        hFigMap.addPlot('hPatch', @patch, XData, YData, vpp, 'EdgeColor', 'k', 'FaceColor', 'flat');
        hFigMap.axApply('default', @axis, [min(XData(:)), max(XData(:)), min(YData(:)), max(YData(:))]);
        hFigMap.axApply('default', @colormap, 'hot');
        hFigMap.addPlot('hText', @text, hCfg.siteLoc(:, 1), hCfg.siteLoc(:, 2), ...
                        arrayfun(@(i) num2str(i), 1:hCfg.nSites, 'UniformOutput', 0), ...
                        'VerticalAlignment', 'bottom', ...
                        'HorizontalAlignment', 'left');
        hFigMap.axApply('default', @xlabel, 'X Position (\mum)');
        hFigMap.axApply('default', @ylabel, 'Y Position (\mum)');
    else
        hFigMap.plotApply('hPatch', @set, 'CData', vpp);
    end

    hFigMap.axApply('default', @title, sprintf('Max: %0.1f \\muVpp', max(clusterVpp)));
    hFigMap.axApply('default', @caxis, [0, max(clusterVpp)]);
end

%% LOCAL FUNCTIONS
function [XData, YData] = getPadCoordinates(hCfg)
    xOffsets = [0 0 1 1] * hCfg.probePad(2);
    yOffsets = [0 1 1 0] * hCfg.probePad(1);
    XData = bsxfun(@plus, hCfg.siteLoc(:, 1)', xOffsets(:));
    YData = bsxfun(@plus, hCfg.siteLoc(:, 2)', yOffsets(:));
end
