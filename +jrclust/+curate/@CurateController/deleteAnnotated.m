function deleteAnnotated(obj)
    %DELETEANNOTATED Delete clusters which are annotated to_delete
    deleteMe = find(strcmp(obj.hClust.clusterNotes, 'to_delete'));
    if ~isempty(deleteMe)
        obj.deleteClusters(deleteMe);
    else
        jrclust.utils.qMsgBox('No units deleted');
    end
end