function uInfo = exportUnitInfo(obj, iCluster)
    %EXPORTUNITINFO Get all data pertinent to a cluster
    uInfo = exportUnitInfo@jrclust.interfaces.Clustering(obj, iCluster);

    if ~isempty(obj.unitSNR)
        uInfo.SNR = obj.unitSNR(iCluster);
    end
end

