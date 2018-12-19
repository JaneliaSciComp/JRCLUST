function hFigSim = doPlotFigSim(hFigSim, hClust, hCfg)
    %DOPLOTFIGSIM Display cluster similarity scores
    hFigSim.figSet('Pointer', 'watch');

    nClusters = hClust.nClusters;
    if isempty(hFigSim.figData) % create from scratch
        hFigSim.axes();
        hFigSim.axSet('Position', [.1 .1 .8 .8], ...
                      'XLimMode', 'manual', ...
                      'YLimMode', 'manual', ...
                      'Layer', 'top', ...
                      'XTick', 1:nClusters, ...
                      'YTick', 1:nClusters);

        if hCfg.fImportKsort
            hFigSim.figData.figView = 'kilosort';
        else
            hFigSim.figData.figView = 'waveform';
        end

        hFigSim.axis([0 nClusters 0 nClusters] + .5);
        hFigSim.axis('xy')
        hFigSim.grid('on');
        hFigSim.xlabel('Cluster #');
        hFigSim.ylabel('Cluster #');

        if strcmp(hFigSim.figData.figView, 'kilosort') && isprop(hClust, 'kSimScore')
            hFigSim.addImagesc('hImSim', 'CData', hClust.kSimScore, hCfg.corrLim);
            hFig.title('[S]plit; [M]erge; [D]elete; [K]iloSort sim score; [W]aveform corr');
        else
            hFigSim.addImagesc('hImSim', 'CData', hClust.simScore, hCfg.corrLim);
            hFigSim.title('[S]plit; [M]erge; [D]elete');
        end

        % selected cluster pair cursors
        hFigSim.addLine('hCursorV', [1 1], [.5 nClusters + .5], 'Color', hCfg.mrColor_proj(2, :), 'LineWidth', 1.5);
        hFigSim.addLine('hCursorH', [.5 nClusters + .5], [1 1], 'Color', hCfg.mrColor_proj(3, :), 'LineWidth', 1.5);

        hFigSim.colorbar();
        hFigSim.addDiag('hDiag', [0, nClusters, 0.5], 'Color', [0 0 0], 'LineWidth', 1.5);
    else
        if strcmp(hFigSim.figData.figView, 'kilosort') && isprop(hClust, 'kSimScore')
            hFigSim.updateImagesc('hImSim', hClust.kSimScore);
            hFigSim.figSet('Name', ['KiloSort cluster similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w')
        else
            hFigSim.updateImagesc('hImSim', hClust.simScore);
            hFigSim.figSet('Name', ['Waveform-based similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w')
        end

        hFigSim.axSet({'XTick', 'YTick'}, {1:nClusters, 1:nClusters});
        hFigSim.addDiag('hDiag', [0, nClusters, 0.5], 'Color', [0 0 0], 'LineWidth', 1.5); % overwrites previous diag plot
    end

    hFigSim.figSet('Pointer', 'arrow');
end
