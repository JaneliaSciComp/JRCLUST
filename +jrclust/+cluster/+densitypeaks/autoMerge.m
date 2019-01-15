function autoMerge(hClust, hCfg, doAssign)
    %AUTOMERGE Automatically merge clusters
    if nargin < 4
        doAssign = 1;
    end

%     if doAssign
%         S_clu = postCluster_(S_clu, P);
%     end

    hClust.refresh();
    hClust.orderClusters('clusterSites');
    hClust.clearNotes();

    if hCfg.outlierThresh > 0
        hClust = rmOutlierSpikes(hClust);
    end

    hClust.doWaveformMerge(hCfg);
    hClust.refresh();
    hClust.orderClusters('clusterSites');
    hClust.updateWaveforms(hCfg);

    hClust.computeCentroids();
    hClust.clearNotes();
    hClust.computeQualityScores();
    commitMsg = sprintf('%s;autoMerge', datestr(now, 31));
    hClust.commit(commitMsg);
end

