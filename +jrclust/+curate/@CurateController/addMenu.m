function addMenu(obj, hFig)
    drawnow;

    outerPosition = hFig.outerPosition;
    hFig.figApply(@set, 'MenuBar','None');

    obj.hMenus('FileMenu') = hFig.uimenu('Label', 'File');
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Save', 'Callback', @(hO, hE) obj.saveFiles());
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Save figures as .fig', 'Callback', @(hO, hE) obj.saveFigures('.fig'));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Save figures as .png', 'Callback', @(hO, hE) obj.saveFigures('.png'));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export units to csv', 'Callback', @(hO, hE) obj.hClust.exportToCSV(), 'Separator', 'on');
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export unit qualities to csv', 'Callback', @(hO, hE) obj.hClust.exportQualityScores());
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export all mean unit waveforms to workspace', 'Callback', @(hO, hE) obj.exportMeanWf(1));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export selected mean unit waveforms to workspace', 'Callback', @(hO, hE) obj.exportMeanWf(0));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export all traces from the selected unit', 'Callback', @(hO, hE) obj.exportTraces());
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export firing rate for all units', 'Callback', @(hO, hE) obj.exportFiringRate());
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Print summary', 'Callback', @(hO, hE) obj.summarize(), 'Separator', 'on');
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Exit', 'Callback', @(hO, hE) obj.endSession(), 'Separator', 'on', 'Accelerator', 'Q');

    obj.hMenus('EditMenu') = hFig.uimenu('Label', 'Edit');
    uimenu(obj.hMenus('EditMenu'), 'Label', '[M]erge', 'Callback', @(hO, hE) obj.mergeSelected());
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Merge auto', 'Callback', @(hO, hE) obj.autoMerge());
    uimenu(obj.hMenus('EditMenu'), 'Label', '[D]elete', 'Callback', @(hO, hE) obj.deleteClusters(), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Delete auto', 'Callback', @(hO, hE) obj.autoDelete());
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Delete annotated', 'Callback', @(hO, hE) obj.deleteAnnotated()); % TW
    uimenu(obj.hMenus('EditMenu'), 'Label', '[S]plit', 'Callback', @(hO, hE) obj.autoSplit(1), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Auto split max-chan', 'Callback', @(hO, hE) obj.autoSplit(0));
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Auto split multi-chan', 'Callback', @(hO, hE) obj.autoSplit(1));
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Recompute selected cluster mean waveform', 'Callback', @(hO, hE) obj.recomputeWaveforms(), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Recompute all mean waveforms', 'Callback', @(hO, hE) obj.recomputeWaveforms(1));
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Reorder clusters by center site', 'Callback', @(hO, hE) obj.reorderClusters('clusterSites'), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Reorder clusters by centroid', 'Callback', @(hO, hE) obj.reorderClusters('Y + X'));

    obj.hMenus('ViewMenu') = uimenu(hFig, 'Label', 'View');
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Show traces', 'Callback', @(hO, hE) obj.showTraces());
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'View all [R]', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'r')));
    uimenu(obj.hMenus('ViewMenu'), 'Label', '[Z]oom selected', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'z')), 'Separator', 'on');
    uimenu(obj.hMenus('ViewMenu'), 'Label', '[W]aveform (toggle)', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'w')));
    uimenu(obj.hMenus('ViewMenu'), 'Label', '[N]umbers (toggle)', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'n')));
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Show raw waveform', 'Callback', @(hO, hE) obj.toggleRaw(hO), 'Separator', 'on');
    %uimenu(obj.hMenus('ViewMenu'), 'Label', 'Threshold by sites', 'Callback', @(hO, hE) keyPressFcn_thresh_(hFig, 'n'));
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Reset window positions', 'Callback', @(hO, hE) obj.resetPositions(), 'Separator', 'on');
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Show config file', 'Callback', @(hO, hE) obj.hCfg.edit(), 'Separator', 'on');

    obj.hMenus('ProjMenu') = uimenu(hFig, 'Label', 'Projection');
    uimenu(obj.hMenus('ProjMenu'), 'Label', 'vpp', 'Callback', @(hO, hE) obj.updateProjection('vpp'));
    uimenu(obj.hMenus('ProjMenu'), 'Label', 'pca', 'Callback', @(hO, hE) obj.updateProjection('pca'));
    uimenu(obj.hMenus('ProjMenu'), 'Label', 'ppca', 'Callback', @(hO, hE) obj.updateProjection('ppca'));
    % uimenu(obj.hMenus('ProjMenu'), 'Label', 'cov', 'Callback', @(hO, hE) obj.updateProjection('cov'));

    obj.hMenus('PlotMenu') = uimenu(hFig, 'Label', 'Plot');
    uimenu(obj.hMenus('PlotMenu'), 'Label', 'Firing rate vs. aux. input (all units)', 'Callback', @(hO, hE) obj.plotAuxRate(0));
    uimenu(obj.hMenus('PlotMenu'), 'Label', 'Firing rate vs. aux. input (selected unit)', 'Callback', @(hO, hE) obj.plotAuxRate(1));

    obj.hMenus('InfoMenu') = uimenu(hFig, 'Label', '', 'Tag', 'InfoMenu');
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Annotate unit', 'Callback', @(hO, hE) obj.annotateUnit('', 1));
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Single unit', 'Callback', @(hO, hE) obj.annotateUnit('single', 0), 'Accelerator', '1');
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Multi unit', 'Callback', @(hO, hE) obj.annotateUnit('multi', 0), 'Accelerator', '2');
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Noise', 'Callback', @(hO, hE) obj.annotateUnit('noise', 0));
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Clear annotation', 'Callback', @(hO, hE) obj.annotateUnit('', 0));
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Equal to', 'Callback', @(hO, hE) obj.annotateUnit('=', 1));

    obj.hMenus('HistMenu') = uimenu(hFig, 'Label', 'History', 'Tag', 'HistMenu');

    obj.hMenus('HelpMenu') = uimenu(hFig, 'Label', 'Help');
    uimenu(obj.hMenus('HelpMenu'), 'Label', '[H]elp', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'h')));

    drawnow;
    hFig.outerPosition = outerPosition;
end