function hFigWav = setFigWavXTicks(hFigWav, hClust, displayCount)
    %SETFIGWAVXTICKS Set X axis ticks for main view
    if displayCount
        xTickLabels = arrayfun(@(i) sprintf('%d (%d)', i, hClust.unitCount(i)), 1:hClust.nClusters, 'UniformOutput', 0);
    else
        xTickLabels = arrayfun(@(i) sprintf('%d', i), 1:hClust.nClusters, 'UniformOutput', 0);
    end

    hFigWav.axApply('default', @set, 'Xtick', 1:hClust.nClusters, ...
                    'XTickLabel', xTickLabels, ...
                    'FontSize', 8);
    if displayCount
        hFigWav.axApply('default', @set, 'XTickLabelRotation', -20);
    else
        hFigWav.axApply('default', @set, 'XTickLabelRotation', 0);
    end

    hFigWav.figData.displayCount = displayCount;
end
