function success = autoMerge(obj)
    %AUTOMERGE Automatically merge clusters
    success = 0;

    if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
        warning('cannot branch from history; use revert() first');
        return;
    end

%     if nargin > 1
%         try
%             obj.hCfg.setTemporaryParams('maxUnitSim', maxUnitSim);
%         catch ME
%             warning('autoMerge aborted: %s', ME.message);
%             return;
%         end
%     end

    obj.refresh(1, []);
    obj.orderClusters('clusterSites');
    obj.clearNotes();

    obj.rmOutlierSpikes();

    obj.updateWaveforms();
    mergedAndUpdated = 0;
    for iRepeat = 1:obj.hCfg.nPassesMerge % single-pass vs dual-pass correction
        nMerged = obj.mergeBySim();
        % recompute mean waveforms and similarity scores to catch what slips
        % through
        if iRepeat > 1 && nMerged < 1 && ~mergedAndUpdated
            obj.updateWaveforms();
            mergedAndUpdated = 1;
        elseif nMerged < 1
            break;
        end
    end

    % no changes, no need to continue
    if iRepeat == 1 && nMerged == 0 && ~isempty(obj.clusterCentroids)
        success = 1;
        return;
    end

    obj.refresh(1, []);
    obj.orderClusters('clusterSites');

    obj.computeCentroids();
    obj.clearNotes();
    obj.computeQualityScores([]);
    obj.commit(sprintf('%s;autoMerge (maxUnitSim %0.2f)', datestr(now, 31), obj.hCfg.maxUnitSim));

    success = 1;
end