function inconsistencies = inconsistentFields(obj)
    %INCONSISTENTFIELDS Check all fields have the correct sizes
    inconsistencies = inconsistentFields@jrclust.interfaces.Clustering(obj);

    % also check all clusterCenters are unique
    if numel(obj.clusterCenters) ~= numel(unique(obj.clusterCenters))
        inconsistencies{end+1} = 'clusterCenters';
    end
end
