function deleteClusters(obj, deleteMe, commitMsg)
    %DELETECLUSTERS Delete clusters either specified or selected
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    if nargin < 2 && numel(obj.selected) > 1
        return;
    elseif nargin < 2
        deleteMe = obj.selected(1);
    end
    if nargin < 3 || isempty(commitMsg)
        deleted = strjoin(arrayfun(@num2str, deleteMe, 'UniformOutput', 0), ', ');
        commitMsg = sprintf('%s;delete;%s', datestr(now, 31), deleted);
    end

    obj.isWorking = 1;
    try
        success = obj.hClust.deleteClusters(deleteMe);
        if success
            % save the new clustering
            obj.hClust.commit(commitMsg);

            obj.isWorking = 0; % in case updateSelect needs to zoom

            if obj.selected > obj.hClust.nClusters
                obj.selected = obj.hClust.nClusters;
            end

            % replot
            obj.updateFigWav();
            obj.updateFigRD(); % centers changed, need replotting
            obj.updateFigSim();
            if numel(deleteMe) == 1 && deleteMe == obj.selected(1)
                obj.updateSelect(deleteMe);
            else
                obj.updateSelect(obj.selected);
            end
        else
            jrclust.utils.qMsgBox('Operation failed.');
        end
    catch ME
        warning('Failed to delete: %s', ME.message)
        jrclust.utils.qMsgBox('Operation failed.');
    end

    obj.isWorking = 0;
end