function annotateUnit(obj, note, doConfirm)
    %ANNOTATEUNIT Add a note to a cluster
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    iCluster = obj.selected(1);

    if nargin < 2
        note = obj.hClust.clusterNotes{iCluster};
    elseif isempty(note) % clear annotation
        obj.hClust.addNote(iCluster, '');
    end

    % set equal to another cluster?
    if ~isempty(note) && strcmp(note, '=') && numel(obj.selected) == 2
        note = sprintf('=%d', obj.selected(2));
    elseif ~isempty(note) && strcmp(note, '=')
        msgbox('Right-click another unit to set equal to selected unit');
        return;
    end

    obj.isWorking = 1;
    try
        if doConfirm
            newNote = inputdlg(sprintf('Cluster %d', iCluster), 'Annotation', 1, {note});
            if ~isempty(newNote)
                obj.hClust.addNote(iCluster, newNote{1});
            end
        else
            obj.hClust.addNote(iCluster, note);
        end

        obj.updateNoteMenu();
    catch ME
        warning('Failed to annotate: %s', ME.message);
        jrclust.utils.qMsgBox('Operation failed.');
    end

    obj.isWorking = 0;
end