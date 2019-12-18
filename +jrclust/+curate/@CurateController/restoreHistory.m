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
                success = 0;

                try
                    success = obj.hClust.revert(entryIndex);
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
                end

            otherwise
                return;
        end
    end
end