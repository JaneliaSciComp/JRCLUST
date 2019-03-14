function uInfo = exportUnitInfo(obj, iCluster)
    %EXPORTUNITINFO Get all data pertinent to a cluster
    if ~ismember(iCluster, 1:obj.nClusters)
        uInfo = [];
        return;
    end

    iSite = obj.clusterSites(iCluster);
    iNeighbors = obj.hCfg.siteNeighbors(:, iSite);

    pos = sprintf('Unit %d (x,y):(%0.1f, %0.1f)[pix]', iCluster, obj.clusterCentroids/obj.hCfg.umPerPix);

    % subsample some (raw or filtered) waveforms
    iSubset = jrclust.utils.subsample(obj.getCenteredSpikes(iCluster), obj.hCfg.nSpikesFigWav);
    if obj.hCfg.showRaw
        meanWf = obj.meanWfGlobalRaw(:, iNeighbors, iCluster);
        if isempty(obj.spikesRawVolt)
            obj.spikesRawVolt = jrclust.utils.rawTouV(obj.spikesRaw, obj.hCfg);
        end
        sampleWf = obj.spikesRawVolt(:, :, iSubset);
        sampleWf = jrclust.filters.fftLowpass(sampleWf, obj.hCfg.getOr('fc_spkwav_show', []), obj.hCfg.sampleRate);
    else
        meanWf = obj.meanWfGlobal(:,iNeighbors,iCluster);
        if isempty(obj.spikesFiltVolt)
            obj.spikesFiltVolt = jrclust.utils.filtTouV(obj.spikesFilt, obj.hCfg);
        end

        sampleWf = obj.spikesFiltVolt(:, :, iSubset);
    end

    uInfo = struct('cluster', iCluster, ...
                   'xyPos', obj.clusterCentroids(iCluster), ...
                   'meanWf', meanWf, ...
                   'neighbors', iNeighbors, ...
                   'position', pos, ...
                   'sampleWf', sampleWf);

    if ~isempty(obj.unitLRatio)
        uInfo.LRatio = obj.unitLRatio(iCluster);
    end
    if ~isempty(obj.unitISIRatio)
        uInfo.ISIRatio = obj.unitISIRatio(iCluster);
    end
    if ~isempty(obj.unitIsoDist)
        uInfo.IsoDist = obj.unitIsoDist(iCluster);
    end
    if ~isempty(obj.unitPeaksRaw)
        uInfo.peaksRaw = obj.unitPeaksRaw(iCluster);
    end
    if ~isempty(obj.unitVppRaw)
        uInfo.vpp = obj.unitVppRaw(iCluster);
    end
end