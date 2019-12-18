function hFigWav = plotFigWav(hFigWav, hClust, maxAmp, showSubset)
    %PLOTFIGWAV Plot the main view (FigWav)
    if nargin < 4
        if ~isfield(hFigWav.figData, 'showSubset')
            hFigWav.figData.showSubset = 1:hClust.nClusters;
        end

        showSubset = hFigWav.figData.showSubset;
    end

    if ~hFigWav.hasAxes('default') || (isfield(hFigWav.figData, 'showSubset') && ~jrclust.utils.isEqual(showSubset, hFigWav.figData.showSubset)) % construct from scratc+h
        hFigWav.addAxes('default');
        hFigWav.axApply('default', @set, 'Position', [.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual');
        hFigWav.axApply('default', @xlabel, 'Unit #');
        hFigWav.axApply('default', @ylabel, 'Site #');
        hFigWav.axApply('default', @grid, 'on');
        hFigWav.axApply('default', @axis, [0, numel(showSubset) + 1, 0, hClust.hCfg.nSites + 1]);

        hFigWav.figData.showSubset = showSubset;

        hFigWav = plotSpikeWaveforms(hFigWav, hClust, maxAmp);
        hFigWav = plotMeanWaveforms(hFigWav, hClust, maxAmp);
        hFigWav.setHideOnDrag('hSpkAll');
    else
        hFigWav = plotSpikeWaveforms(hFigWav, hClust, maxAmp);
        % clear mean waveforms
        iGroup = 1;
        while hFigWav.hasPlot(sprintf('Group%d', iGroup))
            hFigWav.rmPlot(sprintf('Group%d', iGroup));
            iGroup = iGroup + 1;
        end
        % replot mean waveforms
        hFigWav = plotMeanWaveforms(hFigWav, hClust, maxAmp);
    end

    hFigWav.figData.zoom = 1;

    info_ = jrclust.utils.info();
    hFigWav.axApply('default', @title, sprintf('%s v%s; press [H] for help (scale: %0.1f uV)', info_.program, jrclust.utils.version(), maxAmp));
end
