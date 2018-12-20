function hFigWav = setFigWavXTicks(hFigWav, hClust, displayCount)
    %SETFIGWAVXTICKS Set X axis ticks for main view
    if displayCount
        xTickLabels = arrayfun(@(i) sprintf('%d (%d)', i, hClust.clusterCounts(i)), 1:hClust.nClusters, 'UniformOutput', false);
    else
        xTickLabels = arrayfun(@(i) sprintf('%d', i), 1:hClust.nClusters, 'UniformOutput', false);
    end

    hFigWav.axSet('Xtick', 1:hClust.nClusters, ...
                  'XTickLabel', xTickLabels, ...
                  'FontSize', 8);
    if displayCount
        hFigWav.axSet('XTickLabelRotation', -20);
    else
        hFigWav.axSet('XTickLabelRotation', 0);
    end

    hFigWav.figData.displayCount = displayCount;
end