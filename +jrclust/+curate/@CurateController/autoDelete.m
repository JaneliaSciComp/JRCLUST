function autoDelete(obj)
    %AUTODELETE Automatically delete clusters by SNR/spike count
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    hFigDelete = jrclust.views.Figure('', [.5 .7 .35 .3], ['Delete Auto: ', obj.hCfg.sessionName], 0, 0);
    if isfield(obj.hClust,'unitSNR')
        hFigDelete.addPlot('hPlotSNR', obj.hClust.unitSNR(:), obj.hClust.unitCount(:), '.'); % show cluster SNR and spike count
        hFigDelete.axApply('default', @xlabel, 'Unit SNR');
    else
        hFigDelete.addPlot('hPlotSNR', obj.hClust.unitVpp(:), obj.hClust.unitCount(:), '.'); % show cluster SNR and spike count
        hFigDelete.axApply('default', @xlabel, 'Unit \mu Vpp');
    end
    hFigDelete.axApply('default', @ylabel, '# spikes/unit');
    hFigDelete.axApply('default', @grid, 'on');
    hFigDelete.axApply('default', @set, 'YScale', 'log');

    % ask user which clusters to delete
    recDurationSecApprox = obj.hClust.spikeTimes(end)./obj.hCfg.sampleRate;
    % above line uses approximation of number of seconds in the recording based on the time of the last spike.
    % there is a field called "recDurationSec" but it's empty so this is
    % the best I can do. -Adrian Bondy
    unitFR = obj.hClust.unitCount(:)./double(recDurationSecApprox);
    if isfield(obj.hClust,'unitSNR')
        dlgAns = inputdlg({'Min Unit SNR:', 'Max Unit SNR:', sprintf('Min Unit %cVpp:',956),sprintf('Max Unit %cVpp:',956),'Minimum # spikes/unit:','Minimum Firing Rate:'}, 'Auto-deletion based on SNR', 1, {'5', 'inf', '20','inf', '100','0'}); % also ask about # spikes/unit (or firing rate) @TODO
    else
        dlgAns = inputdlg({sprintf('Min Unit %cVpp:',956),sprintf('Max Unit %cVpp:',956),'Minimum # spikes/unit:','Firing Rate:'}, 'Auto-deletion based on SNR', 1, {'20','inf', '100','0'}); % also ask about # spikes/unit (or firing rate) @TODO
    end

    % parse user input
    if isempty(dlgAns)
        closehDelete;
        return;
    end
    if isfield(obj.hClust,'unitSNR')
        snrMin = str2double(dlgAns{1});
        snrMax = str2double(dlgAns{2});
        vppmin = str2double(dlgAns{3});
        vppmax = str2doubl(dlgAns{4});
        minCount = round(str2double(dlgAns{5}));
        minfr = str2double(dlgAns{6});
    else
        vppmin = str2double(dlgAns{1});
        vppmax = str2double(dlgAns{2});
        minCount = round(str2double(dlgAns{3}));
        minfr = str2double(dlgAns{4});
    end

    for i=1:length(dlgAns)
        if isnan(str2double(dlgAns{i}))
            h=jrclust.utils.qMsgBox('Invalid criteria.');
            h.DeleteFcn = @closehDelete;
            return
        end
    end
    if isfield(obj.hClust,'unitSNR')
        deleteMe = find(obj.hClust.unitSNR(:) < snrMin | obj.hClust.unitCount(:) < minCount | obj.hClust.unitSNR(:) > snrMax | unitFR < minfr | obj.hClust.unitVpp < vppmin | obj.hClust.unitVpp > vppmax);
    else
        deleteMe = find( obj.hClust.unitCount(:) < minCount | unitFR < minfr | obj.hClust.unitVpp < vppmin | obj.hClust.unitVpp > vppmax);
    end

    if isempty(deleteMe)
        h=jrclust.utils.qMsgBox('No units deleted.');
        h.DeleteFcn = @closehDelete;
        return;
    end
    if numel(deleteMe) >= obj.hClust.nClusters
        h=jrclust.utils.qMsgBox('Refusing to delete all units.');
        h.DeleteFcn = @closehDelete;
        return;
    end

    % delete and update
    if isfield(obj.hClust,'unitSNR')
        hFigDelete.toForeground;hold on;
        plot(obj.hClust.unitSNR(deleteMe), obj.hClust.unitCount(deleteMe), 'r.'); % show cluster SNR and spike count
    else
        hFigDelete.toForeground;hold on;
        plot( obj.hClust.unitVpp(deleteMe), obj.hClust.unitCount(deleteMe), 'rx'); % show cluster SNR and spike count
    end
    dlgans = questdlg(sprintf('%d units will be deleted. Please confirm.', numel(deleteMe)),'Confirm Deletion','OK','Cancel','OK');
    hFigDelete.close();
    if strcmp(dlgans,'OK')
        obj.deleteClusters(deleteMe);
    else
        return;
    end

    function closehDelete(~,~)
       hFigDelete.close();
    end

end
