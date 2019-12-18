function annotateUnit(obj, note, confirm)
    %ANNOTATEUNIT Add a note to a cluster
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    iCluster = obj.selected(1);

    % set equal to another cluster?
    if strcmp(note, '=') && numel(obj.selected) == 2
        note = sprintf('=%d', obj.selected(2));
    elseif strcmp(note, '=')
        msgbox('Right-click another unit to set equal to selected unit');
        return;
    end

    obj.isWorking = 1;
    try
        if confirm
            newNote = inputdlg(sprintf('Unit %d', iCluster), 'Annotation', 1, {note});
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