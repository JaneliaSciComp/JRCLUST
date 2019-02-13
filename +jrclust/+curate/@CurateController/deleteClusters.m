function deleteClusters(obj, deleteMe)
    %DELETECLUSTERS Delete clusters either specified or selected
    if obj.isWorking
        return;
    end
    obj.isWorking = 1;

    if nargin < 2 && numel(obj.selected) > 1
        return;
    elseif nargin < 2
        deleteMe = obj.selected(1);
    end

    success = obj.hClust.deleteClusters(deleteMe);
    if success
        % save the new clustering
        deleted = strjoin(arrayfun(@num2str, deleteMe, 'UniformOutput', 0), ', ');
        commitMsg = sprintf('%s;delete;%s', datestr(now, 31), deleted);
        obj.hClust.commit(commitMsg);

        % replot
        obj.updateFigWav();
        %obj.updateFigRD(); % centers changed, need replotting
        obj.updateFigSim();
        if numel(deleteMe) == 1 && deleteMe == obj.selected(1)
            obj.updateSelect(deleteMe);
        else
            obj.updateSelect(obj.selected);
        end
    end

    obj.isWorking = 0;
end