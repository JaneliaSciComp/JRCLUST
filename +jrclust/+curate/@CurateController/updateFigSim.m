function updateFigSim(obj)
    %UPDATEFIGSIM Update similarity score figure
    if ~obj.hasFig('FigSim')
        return;
    end

    plotFigSim(obj.hFigs('FigSim'), obj.hClust, obj.hCfg, obj.selected);
end

%% LOCAL FUNCTIONS
function hFigSim = plotFigSim(hFigSim, hClust, hCfg, selected)
    %PLOTFIGSIM Display cluster similarity scores
    hFigSim.wait(1);

    nClusters = hClust.nClusters;
    if ~hFigSim.hasAxes('default') % create from scratch
        hFigSim.addAxes('default');
        hFigSim.axApply('default', @set, 'Position', [.1 .1 .8 .8], ...
                        'XLimMode', 'manual', ...
                        'YLimMode', 'manual', ...
                        'Layer', 'top', ...
                        'XTick', 1:nClusters, ...
                        'YTick', 1:nClusters);

        if isa(hClust, 'jrclust.sort.TemplateClustering')
            hFigSim.figData.figView = 'template'; % start out showing template sim scores
        else
            hFigSim.figData.figView = 'waveform';
        end

        hFigSim.axApply('default', @axis, [0 nClusters 0 nClusters] + .5);
        hFigSim.axApply('default', @axis, 'xy')
        hFigSim.axApply('default', @grid, 'on');
        hFigSim.axApply('default', @xlabel, 'Cluster #');
        hFigSim.axApply('default', @ylabel, 'Cluster #');

        if strcmp(hFigSim.figData.figView, 'template') && isprop(hClust, 'templateSim')
            hFigSim.addPlot('hImSim', @imagesc, 'CData', hClust.templateSim, hCfg.corrRange);
            %hFigSim.axApply('default', @title, '[S]plit; [M]erge; [D]elete; [K]iloSort sim score; [W]aveform corr');
            hFigSim.figApply(@set, 'Name', ['Template-based similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w');
        else
            hFigSim.addPlot('hImSim', @imagesc, 'CData', hClust.waveformSim, hCfg.corrRange);
            %hFigSim.axApply('default', @title, '[S]plit; [M]erge; [D]elete');
        end

        % selected cluster pair cursors
        hFigSim.addPlot('hCursorV', @line, [1 1], [.5 nClusters + .5], 'Color', hCfg.colorMap(2, :), 'LineWidth', 1.5);
        hFigSim.addPlot('hCursorH', @line, [.5 nClusters + .5], [1 1], 'Color', hCfg.colorMap(3, :), 'LineWidth', 1.5);

        hFigSim.axApply('default', @colorbar);
        hFigSim.addDiag('hDiag', [0, nClusters, 0.5], 'Color', [0 0 0], 'LineWidth', 1.5);
    else
        if strcmp(hFigSim.figData.figView, 'template') && isprop(hClust, 'templateSim')
            hFigSim.plotApply('hImSim', @set, 'CData', hClust.templateSim);
            hFigSim.figApply(@set, 'Name', ['Template-based similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w')
        else
            hFigSim.plotApply('hImSim', @set, 'CData', hClust.waveformSim);
            hFigSim.figApply(@set, 'Name', ['Waveform-based similarity score (click): ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w')
        end

        hFigSim.axApply('default', @set, {'XTick', 'YTick'}, {1:nClusters, 1:nClusters});
        hFigSim.addDiag('hDiag', [0, nClusters, 0.5], 'Color', [0 0 0], 'LineWidth', 1.5); % overwrites previous diag plot
    end

    iCluster = selected(1);
    if numel(selected) > 1
        jCluster = selected(2);
    else
        jCluster = iCluster;
    end
    
    if strcmp(hFigSim.figData.figView, 'template')
        scoreij = hClust.templateSim(iCluster, jCluster);
    elseif strcmp(hFigSim.figData.figView, 'waveform')
        scoreij = hClust.waveformSim(iCluster, jCluster);
    end
    hFigSim.axApply('default', @title, sprintf('Cluster %d vs. Cluster %d: %0.3f', iCluster, jCluster, scoreij));

    hFigSim.wait(0);
end