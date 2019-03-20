function success = tryImport(obj, sRes)
    reqFields = {'spikeClusters', 'spikeRho', 'spikeDelta', ...
                 'spikeNeigh', 'ordRho', 'clusterCenters', ...
                 'initialClustering', 'clusterNotes'};
    if all(ismember(reqFields, fieldnames(sRes)))
        obj.spikeClusters = double(sRes.spikeClusters);
        obj.clusterCenters = sRes.clusterCenters;

        
        % object takes initalClustering from sRes.spikeClusters
        if isfield(sRes, 'initialClustering')
            sRes.spikeClusters = double(sRes.initialClustering);
        end

        if ~isfield(sRes, 'ordRho')
            [~, sRes.ordRho] = sort(obj.spikeRho, 'descend');
        end

        % mean waveforms
        meanFields = {'meanWfGlobal', 'meanWfGlobalRaw', 'meanWfLocal', ...
                      'meanWfLocalRaw', 'meanWfRawHigh', 'meanWfRawLow'};
        for i = 1:numel(meanFields)
            mf = meanFields{i};
            if isfield(sRes, mf)
                obj.(mf) = sRes.(mf);
            end
        end

        % cluster-wise summaries
        sumFields = {'unitCount', 'clusterSites', 'spikesByCluster'};
        for i = 1:numel(sumFields)
            mf = sumFields{i};
            if isfield(sRes, mf)
                obj.(mf) = sRes.(mf);
            end
        end

        % quality scores
        qualFields = {'nSitesOverThresh', 'siteRMS', 'unitISIRatio', ...
                      'unitIsoDist', 'unitLRatio', 'unitPeaks', ...
                      'unitPeaksRaw', 'unitSNR', 'unitVpp', 'unitVppRaw'};
        for i = 1:numel(qualFields)
            mf = qualFields{i};
            if isfield(sRes, mf)
                obj.(mf) = sRes.(mf);
            end
        end

        % sim score
        if isfield(sRes, 'waveformSim')
            obj.waveformSim = sRes.waveformSim;
        else
            obj.computeWaveformSim([]);
        end

        if isfield(sRes, 'clusterNotes')
            obj.clusterNotes = sRes.clusterNotes;
        else
            obj.clearNotes();
        end

        obj.sRes = sRes;
        % any mean waveform fields empty? recompute
        if any(cellfun(@(f) isempty(obj.(f)), meanFields))
            obj.updateWaveforms([]);
        end

        % any summary fields empty? refresh
        if any(cellfun(@(f) isempty(obj.(f)), sumFields))
            obj.refresh(1, []);
        end

        % any quality scores empty? recompute
        if any(cellfun(@(f) isempty(obj.(f)), qualFields))
            obj.computeQualityScores([]);
        end
        commitMsg = sprintf('%s;initial import', datestr(now, 31));
        obj.commit(commitMsg);

        success = 1;
    else
        success = 0;
    end
end