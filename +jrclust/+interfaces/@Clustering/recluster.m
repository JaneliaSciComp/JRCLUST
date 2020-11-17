function success = recluster(obj)
%RECLUSTER Reassign clusters, e.g., after a change of parameters.
success = 0;

try
    obj.sRes = jrclust.sort.assignClusters(obj.dRes, obj.sRes, obj.hCfg);
    obj.doRecompute();
    [success, ~] = obj.reorderBy('clusterSites');
    success = 1;
catch ME
    warning('Failed to reassign: %s', ME.message);
end
end %fun
