function autoDelete(obj)
    %AUTODELETE Automatically delete clusters by SNR/spike count
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    elseif ~isa(obj.hClust, 'jrclust.sort.DensityPeakClustering')
        jrclust.utils.qMsgBox('Operation not supported for this type of clustering');
        return;
    end

    hFigDelete = jrclust.views.Figure('', [.5 .7 .35 .3], ['Delete Auto: ', obj.hCfg.sessionName], 0, 0);

    hFigDelete.addPlot('hPlotSNR', obj.hClust.unitSNR(:), obj.hClust.unitCount(:), '.'); % show cluster SNR and spike count
    hFigDelete.axApply('default', @xlabel, 'Unit SNR');
    hFigDelete.axApply('default', @ylabel, '# spikes/unit');
    hFigDelete.axApply('default', @grid, 'on');
    hFigDelete.axApply('default', @set, 'YScale', 'log');

    % ask user which clusters to delete
    dlgAns = inputdlg({'Min Unit SNR:', 'Max Unit SNR:', 'Minimum # spikes/unit'}, 'Auto-deletion based on SNR', 1, {'5', 'inf', '0'}); % also ask about # spikes/unit (or firing rate) @TODO
    hFigDelete.close();

    % parse user input
    if isempty(dlgAns)
        return;
    end

    snrMin = str2double(dlgAns{1});
    snrMax = str2double(dlgAns{2});
    minCount = round(str2double(dlgAns{3}));

    if any(isnan([snrMin, snrMax, minCount]))
        jrclust.utils.qMsgBox('Invalid criteria.');
        return;
    end

    deleteMe = find(obj.hClust.unitSNR(:) < snrMin | obj.hClust.unitCount(:) < minCount | obj.hClust.unitSNR(:) > snrMax);
    if isempty(deleteMe)
        jrclust.utils.qMsgBox('No units deleted.');
        return;
    end
    if numel(deleteMe) >= obj.hClust.nClusters
        jrclust.utils.qMsgBox('Cannot delete all units.');
        return;
    end

    % delete and update
    commitMsg = sprintf('%s;autodelete %d (min SNR %0.1f, max SNR %0.1f, %d spikes/unit minimum);%s', datestr(now, 31), numel(deleteMe), snrMin, snrMax, minCount, deleteMe);
    obj.deleteClusters(deleteMe, commitMsg);
    jrclust.utils.qMsgBox(sprintf('Deleted %d units (min SNR %0.1f, max SNR %0.1f, %d spikes/unit minimum).', numel(deleteMe), snrMin, snrMax, minCount));
end