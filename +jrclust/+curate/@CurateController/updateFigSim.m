function updateFigSim(obj)
    %UPDATEFIGSIM Update similarity score figure
    if ~obj.hasFig('FigSim')
        return;
    end

    plotFigSim(obj.hFigs('FigSim'), obj.hClust, obj.hCfg, obj.selected, obj.showSubset);
end

%% LOCAL FUNCTIONS
function hFigSim = plotFigSim(hFigSim, hClust, hCfg, selected, showSubset)
    %PLOTFIGSIM Display cluster similarity scores
    if ~isfield(hFigSim.figData, 'showSubset')
        hFigSim.figData.showSubset = showSubset;
    end

    hFigSim.wait(1);

    xyLabels = arrayfun(@num2str, showSubset, 'UniformOutput', 0);
    nClusters = numel(showSubset);
    if ~hFigSim.hasAxes('default') || ~jrclust.utils.isEqual(showSubset, hFigSim.figData.showSubset) % create from scratch
        hFigSim.addAxes('default');
        hFigSim.axApply('default', @set, 'Position', [.1 .1 .8 .8], ...
                        'XLimMode', 'manual', ...
                        'YLimMode', 'manual', ...
                        'Layer', 'top', ...
                        'XTick', 1:nClusters, ...
                        'XTickLabel', xyLabels, ...
                        'YTick', 1:nClusters, ...
                        'YTickLabel', xyLabels);
        
        if ~isfield(hFigSim.figData,'figView')
            if isa(hClust, 'jrclust.sort.TemplateClustering')
                hFigSim.figData.figView = 'template'; % start out showing template sim scores
            else
                hFigSim.figData.figView = 'waveform';
            end
        end

        hFigSim.axApply('default', @axis, [0 nClusters 0 nClusters] + .5);
        hFigSim.axApply('default', @axis, 'xy')
        hFigSim.axApply('default', @grid, 'on');
        hFigSim.axApply('default', @xlabel, 'Unit #');
        hFigSim.axApply('default', @ylabel, 'Unit #');

        if strcmp(hFigSim.figData.figView, 'template') && isprop(hClust, 'templateSim')
            hFigSim.addPlot('hImSim', @imagesc, 'CData', hClust.templateSim(showSubset, showSubset), hCfg.corrRange);
            hFigSim.figApply(@set, 'Name', ['Template-based similarity score: ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w');
        else
            hFigSim.addPlot('hImSim', @imagesc, 'CData', hClust.waveformSim(showSubset, showSubset), hCfg.corrRange);
        end

        % selected cluster pair cursors
        hFigSim.addPlot('hCursorV', @line, [1 1], [.5 nClusters + .5], 'Color', hCfg.colorMap(2, :), 'LineWidth', 1.5);
        hFigSim.addPlot('hCursorH', @line, [.5 nClusters + .5], [1 1], 'Color', hCfg.colorMap(3, :), 'LineWidth', 1.5);

        hFigSim.axApply('default', @colorbar);
        hFigSim.addDiag('hDiag', [0, nClusters, 0.5], 'Color', [0 0 0], 'LineWidth', 1.5);
    else
        if strcmp(hFigSim.figData.figView, 'template') && isprop(hClust, 'templateSim')
            hFigSim.plotApply('hImSim', @set, 'CData', hClust.templateSim(showSubset, showSubset));
            hFigSim.figApply(@set, 'Name', ['Template-based similarity score: ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w')
        else
            hFigSim.plotApply('hImSim', @set, 'CData', hClust.waveformSim(showSubset, showSubset));
            hFigSim.figApply(@set, 'Name', ['Waveform-based similarity score: ', hCfg.sessionName], 'NumberTitle', 'off', 'Color', 'w')
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
    hFigSim.axApply('default', @title, sprintf('Unit %d vs. Unit %d: %0.3f', iCluster, jCluster, scoreij));

    hFigSim.wait(0);
end