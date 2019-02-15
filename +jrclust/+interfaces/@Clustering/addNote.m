function addNote(obj, iCluster, note)
    %ADDNOTE Annotate a cluster
    if iCluster < 1 || iCluster > obj.nClusters
        return;
    end

    obj.clusterNotes{iCluster} = note;
end