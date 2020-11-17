function success = reassign(obj)
%REASSIGN Reassign clusters, e.g., after a change of parameters
success = 0;

obj.sRes = jrclust.sort.assignClusters(obj.dRes, obj.sRes, obj.hCfg);
obj.commit(obj.spikeClusters, struct(), 'reassign clusters');
obj.orderClusters('clusterSites');
success = 1;
end %fun