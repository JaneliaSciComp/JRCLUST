function hFigWav = plotFigWav(hFigWav, hClust, maxAmp, showSubset, channel_idx)
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

        % hFigWav.figData.vcTitle = 'Scale: %0.1f uV; [H]elp; [Left/Right]:Select cluster; (Sft)[Up/Down]:scale; [M]erge; [S]plit auto; [D]elete; [A]:Resample spikes; [P]STH; [Z]oom; in[F]o; [Space]:Find similar [0]:Annotate Delete [1]:Annotate Signle [2]:Annotate Multi'; % TW
        % hFigWav.axApply('default', @title, sprintf('Scale: %0.1f uV; [H]elp; [Left/Right]:Select cluster; (Sft)[Up/Down]:scale; [M]erge; [S]plit auto; [D]elete; [A]:Resample spikes; [P]STH; [Z]oom; in[F]o; [Space]:Find similar [0]:Annotate Delete [1]:Annotate Signle [2]:Annotate Multi', maxAmp));

        hFigWav.axApply('default', @axis, [0, hClust.nClusters + 1, 0, hClust.hCfg.nSites + 1]);
        hFigWav = plotSpikeWaveforms(hFigWav, hClust, maxAmp, channel_idx);
        hFigWav = plotMeanWaveforms(hFigWav, hClust, maxAmp, channel_idx);
        hFigWav.setHideOnDrag('hSpkAll');
    else
        hFigWav = plotSpikeWaveforms(hFigWav, hClust, maxAmp, channel_idx);
        % clear mean waveforms
        iGroup = 1;
        while hFigWav.hasPlot(sprintf('Group%d', iGroup))
            hFigWav.rmPlot(sprintf('Group%d', iGroup));
            iGroup = iGroup + 1;
        end
        % replot mean waveforms
        hFigWav = plotMeanWaveforms(hFigWav, hClust, maxAmp,channel_idx);
    end

    hFigWav.figData.zoom = 1;

    info_ = jrclust.utils.info();
    hFigWav.axApply('default', @title, sprintf('%s v%s; press [H] for help (scale: %0.1f uV)', info_.program, jrclust.utils.version(), maxAmp));
end
