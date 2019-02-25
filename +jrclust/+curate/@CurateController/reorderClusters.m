function reorderClusters(obj, by)
    %REORDERCLUSTERS Reorder clusters by centroids
    if nargin < 2
        by = 'clusterSites';
    end

    argsort = obj.hClust.orderClusters(by);

    if issorted(argsort)
        jrclust.utils.qMsgBox('Clusters already in order');
    else
        nChanged = sum(argsort(:) ~= (1:obj.hClust.nClusters)');
        jrclust.utils.qMsgBox(sprintf('%d clusters changed', nChanged));
    end

    obj.updateFigWav();
    obj.updateFigSim();
    obj.updateSelect(1); 
end

