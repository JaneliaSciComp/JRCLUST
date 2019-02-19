function restoreHistory(obj, entryIndex)
    %RESTOREHISTORY Restore a position in hClust history
    if entryIndex <= obj.hClust.nEdits
        if isprop(obj.hCfg, 'testRun') && obj.hCfg.testRun
            dlgAns = 'Yes';
        else
            qstring = 'This operation is irreversible. Continue?';
            dlgAns = questdlg(qstring, 'Revert history', 'No');
        end

        switch dlgAns
            case 'Yes'
                hBox = jrclust.utils.qMsgBox('Reverting your history... (this closes automatically)', 0, 1);
                obj.isWorking = 1;
                try
                    obj.hClust.revert(entryIndex-1);
                catch ME
                    warning('Failed to revert: %s', ME.message);
                end
                obj.isWorking = 0; % in case updateSelect needs to zoom

                % success; replot
                jrclust.utils.tryClose(hBox);
                obj.updateFigWav();
                obj.updateFigRD(); % centers changed, need replotting
                obj.updateFigSim();
                obj.updateSelect(1);

            otherwise
                return;
        end
    end
end