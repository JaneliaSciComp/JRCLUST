function deleteAnnotated(obj)
    %DELETEANNOTATED Delete clusters which are annotated to_delete
    deleteMe = find(strcmp(obj.hClust.clusterNotes, 'to_delete'));
    if ~isempty(deleteMe)
        obj.deleteClusters(deleteMe);
    end
end