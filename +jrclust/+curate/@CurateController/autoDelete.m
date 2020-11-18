function success = autoDelete(obj)
%AUTODELETE Automatically delete clusters by SNR/spike count
success = 0;

if obj.isWorking
    jrclust.utils.qMsgBox('An operation is in progress.');
    return;
end

hFigDelete = jrclust.views.Figure('', [.5 .7 .35 .3], ['Delete Auto: ', obj.hCfg.sessionName], 0, 0);
if isprop(obj.hClust, 'unitSNR')
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
if isprop(obj.hClust, 'unitSNR')
    prompt = {'Min Unit SNR:', 'Max Unit SNR:', sprintf('Min Unit %cVpp:', 956), sprintf('Max Unit %cVpp:', 956), 'Minimum # spikes/unit:', 'Minimum Firing Rate:'};
    definput = {'5', 'inf', '20','inf', '100','0'};
else
    prompt = {sprintf('Min Unit %cVpp:', 956), sprintf('Max Unit %cVpp:', 956), 'Minimum # spikes/unit:', 'Minimum Firing Rate:'};
    definput = {'20', 'inf', '100', '0'};
end

if obj.hCfg.getOr('testRun', 0)
    dlgAns = definput;
else
    dlgAns = inputdlg(prompt, 'Auto-deletion based on SNR', 1, definput);
end

% parse user input
if isempty(dlgAns)
    closehDelete;
    return;
end

if isprop(obj.hClust, 'unitSNR')
    snrMin = str2double(dlgAns{1});
    snrMax = str2double(dlgAns{2});
    vppMin = str2double(dlgAns{3});
    vppMax = str2double(dlgAns{4});
    minCount = round(str2double(dlgAns{5}));
    minFR = str2double(dlgAns{6});
else
    vppMin = str2double(dlgAns{1});
    vppMax = str2double(dlgAns{2});
    minCount = round(str2double(dlgAns{3}));
    minFR = str2double(dlgAns{4});
end

% check all values are numeric
if any(cellfun(@(x) isnan(str2double(x)), dlgAns))
    h = jrclust.utils.qMsgBox('Invalid criteria.');
    h.DeleteFcn = @closehDelete;
    return;
end

% look for low unit counts, low firing rates, or Vpp outside desired range
deleteMask = obj.hClust.unitCount(:) < minCount | unitFR < minFR | ...
    obj.hClust.unitVpp < vppMin | obj.hClust.unitVpp > vppMax;

if isprop(obj.hClust, 'unitSNR')
    % also look for SNR outside desired range
    deleteMask = deleteMask | obj.hClust.unitSNR(:) < snrMin | ...
        obj.hClust.unitSNR(:) > snrMax;
end

deleteMe = find(deleteMask);

if isempty(deleteMe)
    h = jrclust.utils.qMsgBox('No units deleted.');
    h.DeleteFcn = @closehDelete;
    success = 1;
    return;
end

if numel(deleteMe) >= obj.hClust.nClusters
    h = jrclust.utils.qMsgBox('Refusing to delete all units.');
    h.DeleteFcn = @closehDelete;
    return;
end

% show units to be deleted
hFigDelete.toForeground; hold on;
if isprop(obj.hClust, 'unitSNR')
    % show cluster SNR and spike count
    xVal = obj.hClust.unitSNR(deleteMe);
else
    % show cluster Vpp and spike count
    xVal = obj.hClust.unitVpp(deleteMe);
end
plot(xVal, obj.hClust.unitCount(deleteMe), 'rx');

% confirm deletion
if obj.hCfg.getOr('testRun', 0)
    dlgans = 'OK';
else
    dlgans = questdlg(sprintf('%d units will be deleted. Please confirm.', ...
        numel(deleteMe)), 'Confirm Deletion', 'OK', 'Cancel', 'OK');
end
hFigDelete.close();

if strcmp(dlgans, 'OK')
    success = obj.deleteClusters(deleteMe);
else
    success = 1;
    return;
end

function closehDelete(~,~)
   hFigDelete.close();
end
end %fun
