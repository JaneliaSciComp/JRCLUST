function success = autoMerge(obj, maxUnitSim)
%AUTOMERGE Automatically merge units based on their similarity scores.
success = 0;

if obj.isWorking
    jrclust.utils.qMsgBox('An operation is in progress.');
    return;
end

if nargin < 2
    if obj.hCfg.getOr('testRun', 0)
        dlgAns = {num2str(obj.hCfg.maxUnitSim)};
    else
        dlgAns = inputdlg('Waveform correlation threshold (0-1):', 'Auto-merge based on waveform threshold', 1, {num2str(obj.hCfg.maxUnitSim)});
    end
    % parse user input
    if isempty(dlgAns)
        return;
    end

    maxUnitSim = str2double(dlgAns{1});
    if isnan(maxUnitSim) || maxUnitSim <= 0 || maxUnitSim > 1
        jrclust.utils.qMsgBox('Invalid criteria.');
        return;
    end
end

if obj.hasFig('FigWav')
    hFigWav = obj.hFigs('FigWav');
    hFigWav.wait(1);
end

nClustersPrev = obj.hClust.nClusters;
hBox = jrclust.utils.qMsgBox('Merging... (this closes automatically)', 0, 1);

obj.isWorking = 1;
try
    obj.hClust.hCfg.setTemporaryParams('maxUnitSim', maxUnitSim);
    success = obj.hClust.autoMerge();
    if success % replot
        jrclust.utils.tryClose(hBox);

        obj.isWorking = 0; % in case updateSelect needs to zoom

        obj.updateFigWav();
        obj.updateFigRD(); % centers changed, need replotting
        obj.updateFigSim();
        obj.updateSelect(obj.selected, 1);

        jrclust.utils.qMsgBox(sprintf('Merged %d clusters >%0.2f maxUnitSim.', nClustersPrev - obj.hClust.nClusters, maxUnitSim));

        obj.hClust.doRecompute();
    else
        jrclust.utils.tryClose(hBox);
        jrclust.utils.qMsgBox('Auto merge failed.');
    end

    obj.hClust.hCfg.resetTemporaryParams('maxUnitSim');

    if obj.hasFig('FigWav')
        hFigWav = obj.hFigs('FigWav');
        hFigWav.wait(0);
    end
catch ME
    jrclust.utils.tryClose(hBox);
    warning('Failed to merge: %s', ME.message);
    jrclust.utils.qMsgBox('Operation failed.');
    rethrow(ME)
end

obj.isWorking = 0;
end