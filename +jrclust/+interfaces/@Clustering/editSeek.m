function editSeek(obj, seekTo)
    %EDITSEEK Seek back and forth to position in history
    if seekTo < 0 || seekTo > obj.nEdits
        return;
    end

    % restore initial state and fast forward
    spikeClusters_ = obj.initialClustering;
    if isfield(obj.sRes, 'clusterCenters')
        obj.clusterCenters = obj.sRes.clusterCenters;
    else
        obj.clusterCenters = [];
    end

    for j = 1:seekTo+1
        op = obj.history(j, :);
        if strcmp(op{2}, 'delete')
            deleted = op{3};
            spikeClusters_(ismember(spikeClusters_, deleted)) = 0;
        elseif strcmp(op{2}, 'merge')
            iCluster = op{3};
            jCluster = op{4};
            spikeClusters_(spikeClusters_ == jCluster) = iCluster;

            % shift clusters larger than jCluster down by 1
            gtMask = (spikeClusters_ > jCluster);
            spikeClusters_(gtMask) = spikeClusters_(gtMask) - 1;
        elseif strcmp(op{2}, 'split')
            iCluster = op{3};
            retained = op{4};

            % shift clusters larger than iCluster up by 1 (make
            % room for splitted off cluster)
            gtMask = (spikeClusters_ > iCluster);
            spikeClusters_(gtMask) = spikeClusters_(gtMask) + 1;

            % take splitted off spikes and make a new cluster
            % of them
            newMask = (spikeClusters_ == iCluster & ~ismember(spikeClusters_, retained));
            spikeClusters_(newMask) = iCluster + 1;
        else
            diffs = obj.history{j, 3};
            if isempty(diffs)
                continue;
            end

            iDiffs = diffs(1, :);
            sDiffs = diffs(2, :);
            spikeClusters_(iDiffs) = sDiffs;
        end
    end
    obj.spikeClusters = spikeClusters_;

    obj.refresh(0, []);
    if ~isempty(obj.meanWfGlobal) % only compute mean waveforms if we already had them
        obj.updateWaveforms();
    end

    if ~isempty(obj.unitVpp)
        obj.computeQualityScores([]);
    end
end