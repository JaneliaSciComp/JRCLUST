function clearNotes(obj)
    %CLEARNOTES Remove all cluster notes
    obj.clusterNotes = arrayfun(@(~) '', 1:obj.nClusters, 'UniformOutput', 0);
end