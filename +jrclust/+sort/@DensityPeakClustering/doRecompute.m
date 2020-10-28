function success = doRecompute(obj)
%DORECOMPUTE Recompute derivative properties of the spike table after an
%edit.
%   After basic doRecompute, update clusterCenters, then clear
%   obj.recompute.
success = doRecompute@jrclust.interfaces.Clustering(obj);

% recompute clusterCenters
if success
    success = obj.assignClusterCenters();
end

if success
    obj.recompute = [];
end
end % fun
