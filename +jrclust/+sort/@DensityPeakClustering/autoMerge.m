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
    for iRepeat = 1:obj.hCfg.nPassesMerge % single-pass vs dual-pass correction
        nMerged = obj.mergeBySim();
        if nMerged < 1
            break;
        end
    end

    obj.refresh(1, []);
    obj.orderClusters('clusterSites');
    obj.updateWaveforms();

    obj.computeCentroids();
    obj.clearNotes();
    obj.computeQualityScores([]);
    obj.commit(sprintf('%s;autoMerge', datestr(now, 31)));

    if nargin == 2
        obj.hCfg.resetTemporaryParams('maxUnitSim');
    end

    success = 1;
end