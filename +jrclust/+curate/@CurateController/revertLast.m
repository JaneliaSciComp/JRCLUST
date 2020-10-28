function revertLast(obj, nReverts)
%REVERTLAST Revert the last `nReverts` operations in the history table.
if nReverts > obj.hClust.nEdits || nReverts <= 0
    return;
end

if isprop(obj.hCfg, 'testRun') && obj.hCfg.testRun
    dlgAns = 'Yes';
else
    dlgAns = questdlg('Are you sure you wish to revert?', 'Revert history', 'No');
end

switch dlgAns
    case 'Yes'
        hBox = jrclust.utils.qMsgBox('Reverting your history... (this closes automatically)', 0, 1);
        obj.isWorking = 1;
        success = 0;

        try
            success = obj.hClust.revertLast(nReverts);
            obj.showSubset = 1:obj.hClust.nClusters;
        catch ME
            warning('Failed to revert: %s', ME.message);
        end

        obj.isWorking = 0; % in case updateSelect needs to zoom
        jrclust.utils.tryClose(hBox);

        % success; replot
        if success
            obj.updateFigWav();
            obj.updateFigRD(); % centers changed, need replotting
            obj.updateFigSim();
            obj.updateSelect(1);
            obj.updateHistMenu();
            obj.keyPressFigWav([], struct('Key', 'z')); % zoom
        end

    otherwise
        return;
end
end