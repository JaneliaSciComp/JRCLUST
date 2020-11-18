function success = reorderClusters(obj, by)
%REORDERCLUSTERS Reorder clusters by some criterion.
if nargin < 2
    by = 'clusterSites';
end

[success, argsort] = obj.hClust.reorderBy(by);

if ~success
    jrclust.utils.qMsgBox('Operation failed.');
    return
elseif issorted(argsort)
    jrclust.utils.qMsgBox('Clusters already in order');
    return;
else
    nChanged = sum(argsort(:) ~= (1:obj.hClust.nClusters)');
    jrclust.utils.qMsgBox(sprintf('%d clusters changed', nChanged));
    obj.showSubset = find(ismember(argsort, obj.showSubset));
end

% select the currently selected unit at its new location
obj.selected = find(obj.selected == argsort);
obj.updateFigWav();
obj.updateFigSim();
obj.updateSelect(obj.selected, 1); 
end %fun

