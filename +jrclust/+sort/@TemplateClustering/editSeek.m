function editSeek(obj, seekTo)
    %EDITSEEK Seek back and forth to position in history
    if seekTo < 0 || seekTo > obj.nEdits
        return;
    end

    % restore initial state and fast forward
    spikeClusters_ = obj.initialClustering;

    for j = 1:seekTo+1
        op = obj.history(j, :);
        switch op{2}
            case 'delete'
                deleted = op{3};
                nClustersOld = max(spikeClusters_);
                spikeClusters_(ismember(spikeClusters_, deleted)) = 0;
                if ~isempty(obj.clusterCenters)
                    obj.clusterCenters(deleted) = [];
                end

                % shift clusters larger than deleted down by 1
                if numel(deleted) == 1
                    gtMask = (spikeClusters_ > deleted);
                    spikeClusters_(gtMask) = spikeClusters_(gtMask) - 1;
                else
                    keepMe = setdiff(1:nClustersOld, deleted);
                    nClustersNew = numel(keepMe);

                    good = (spikeClusters_ > 0);
                    mapFrom = zeros(1, nClustersOld);
                    mapFrom(keepMe) = 1:nClustersNew;

                    spikeClusters_(good) = mapFrom(spikeClusters_(good));
                end

            case 'merge'
                iCluster = op{3};
                jCluster = op{4};
                spikeClusters_(spikeClusters_ == jCluster) = iCluster;

                % shift clusters larger than jCluster down by 1
                gtMask = (spikeClusters_ > jCluster);
                spikeClusters_(gtMask) = spikeClusters_(gtMask) - 1;

            case 'split'
                iCluster = op{3};
                retained = op{4};

                % shift clusters larger than iCluster up by 1 (make
                % room for splitted off cluster)
                gtMask = (spikeClusters_ > iCluster);
                spikeClusters_(gtMask) = spikeClusters_(gtMask) + 1;

                % take splitted off spikes and make a new cluster
                % of them
                iMask = (spikeClusters_ == iCluster);
                iSpikes = find(iMask);
                jSpikes = iSpikes(~ismember(iSpikes, retained)); % spikes to split off
                spikeClusters_(jSpikes) = iCluster + 1;

            case 'partition'
                iCluster = op{3};
                assignPart = op{4};
                initialSpikes = find(spikeClusters_ == iCluster);

                for kSplit = 1:numel(assignPart)-1
                    if kSplit > 1
                        iSpikes = find(spikeClusters_ == iCluster);
                    else
                        iSpikes = initialSpikes;
                    end

                    splitOff = initialSpikes(assignPart{kSplit});
                    retain = iSpikes(~ismember(iSpikes, splitOff));

                    iSite = mode(obj.spikeSites(retain));
                    jSite = mode(obj.spikeSites(splitOff));

                    % swap retain and split-off spikes; always make splitOff the next
                    % cluster
                    if iSite > jSite
                        splitOff = retain;
                    end

                    % make room for new cluster
                    mask = (spikeClusters_ > iCluster + kSplit - 1);
                    spikeClusters_(mask) = spikeClusters_ + 1;
                    spikeClusters_(splitOff) = iCluster + kSplit;
                end

            otherwise
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