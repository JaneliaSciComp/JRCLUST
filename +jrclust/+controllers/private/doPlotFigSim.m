function hFigSim = doPlotFigSim(hFigSim, hClust, hCfg)
    %DOPLOTFIGSIM Display cluster similarity scores
    hFigSim.figSet('Pointer', 'watch');

    nClusters = hClust.nClusters;
    if isempty(hFigSim.figMetadata)
        % figMetadata.hAx = axes_new_(hFig);
        hFigSim.axes();
        % set(figMetadata.hAx, 'Position', [.1 .1 .8 .8], 'XLimMode', 'manual', 'YLimMode', 'manual', 'Layer', 'top');
        hFigSim.axSet('Position', [.1 .1 .8 .8], 'XLimMode', 'manual', 'YLimMode', 'manual', 'Layer', 'top');
        % set(figMetadata.hAx, {'XTick', 'YTick'}, {1:nClusters, 1:nClusters});
        hFigSim.axSet({'XTick', 'YTick'}, {1:nClusters, 1:nClusters});

        if hCfg.fImportKsort
            hFigSim.figMetadata.figView = 'kilosort';
        else
            hFigSim.figMetadata.figView = 'waveform';
        end

        % axis_(figMetadata.hAx, [0 nClusters 0 nClusters] + .5);
        hFigSim.axis([0 nClusters 0 nClusters] + .5);
        % axis(figMetadata.hAx, 'xy');
        hFigSim.axis('xy')
        % grid(figMetadata.hAx, 'on');
        hFigSim.grid('on');
        %xlabel(figMetadata.hAx, 'Clu#');
        hFigSim.xlabel('Cluster #');
        %ylabel(figMetadata.hAx, 'Clu#');
        hFigSim.ylabel('Cluster #');

        %clears title and current figure
        if strcmp(hFigSim.figMetadata.figView, 'kilosort') && isprop(hClust, 'kSimScore')
            % figMetadata.hImWavCor = imagesc(hClust.mrSim_clu, hCfg.corrLim); 
            hFigSim.figMetadata.hImSim = imagesc(hClust.kSimScore, hCfg.corrLim);
            hFig.title('[S]plit; [M]erge; [D]elete; [K]iloSort sim score; [W]aveform corr');
        else
            % figMetadata.hImWavCor = imagesc(hClust.mrWavCor, hCfg.corrLim); %clears title and current figure
            hFigSim.figMetadata.hImSim = imagesc(hClust.simScore, hCfg.corrLim);
            hFigSim.title('[S]plit; [M]erge; [D]elete');
        end

        hFigSim.figMetadata.hCursorV = line([1 1], [.5 nClusters + .5], 'Color', [0 0 0], 'LineWidth', 1.5);
        hFigSim.figMetadata.hCursorH = line([.5 nClusters + .5], [1 1], 'Color', [1 0 0], 'LineWidth', 1.5);
        % colorbar(figMetadata.hAx);
        hFigSim.colorbar();

%         if get_set_(hCfg, 'fImportKsort', 0)
%             figMetadata.vcTitle = '[S]plit; [M]erge; [D]elete; [K]iloSort sim score; [W]aveform corr';
%         else
%             figMetadata.vcTitle = '[S]plit; [M]erge; [D]elete';
%         end

        % set(hFigSim, 'KeyPressFcn', @keyPressFcn_FigWavCor_);
        % mouse_figure(hFigSim, figMetadata.hAx, @button_FigWavCor_);
        hFigSim.figMetadata.hDiag = plotDiag_([0, nClusters, .5], 'Color', [0 0 0], 'LineWidth', 1.5);
    else
        if strcmp(hFigSim.figMetadata.figView, 'kilosort') && isprop(hClust, 'kSimScore')
            set(hFigSim.figMetadata.hImWavCor, 'CData', hClust.kSimScore);
            % set(hFigSim, 'Name', ['KiloSort cluster similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w');
            hFigSim.figSet('Name', ['KiloSort cluster similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w')
        else
            set(hFigSim.figMetadata.hImWavCor, 'CData', hClust.mrWavCor);
            % set(hFigSim, 'Name', ['Waveform-based similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w');
            hFigSim.figSet('Name', ['Waveform-based similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w')
        end

        %set(figMetadata.hAx, {'XTick', 'YTick'}, {1:nClusters, 1:nClusters});
        hFigSim.axSet({'XTick', 'YTick'}, {1:nClusters, 1:nClusters});
        [diagX, diagY] = plotDiag__([0, nClusters, .5]);
        set(hFigSim.figMetadata.hDiag, 'XData', diagX, 'YData', diagY);
    end

    % output
    % set(hFigSim, 'UserData', figMetadata);
    hFigSim.figSet('Pointer', 'arrow');
end
