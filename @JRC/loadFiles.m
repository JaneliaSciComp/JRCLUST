function loadFiles(obj)
%LOADFILES Load results struct and binary files.
obj.loadRes();
obj.loadBinaries();

if isfield(obj.res, 'spikeTimes')
    if isfield(obj.res, 'history') && iscell(obj.res.history) % old-style history
        obj.res.history = convertHistory(obj.res.history, obj.res.initialClustering, obj.hCfg);
    end
else
    warning('spikeTimes not found in %s', obj.hCfg.resFile);
    res_ = struct();
end

obj.res = res_;
end

function history = convertHistory(oldHistory, spikeClusters, hCfg)
%CONVERTHISTORY Convert old-style (cell) history to the history file
history = containers.Map('KeyType', 'int32', 'ValueType', 'char');
history(1) = 'initial commit';

fidHist = fopen(hCfg.histFile, 'w');
fwrite(fidHist, int32(1), 'int32');
fwrite(fidHist, int32(spikeClusters), 'int32');

for j = 2:size(oldHistory, 1)
    op = oldHistory(j, :);

    switch op{2}
        case 'delete'
            deleted = op{3};
            nClustersOld = max(spikeClusters);
            spikeClusters(ismember(spikeClusters, deleted)) = 0;
            if ~isempty(obj.clusterCenters)
                obj.clusterCenters(deleted) = [];
            end

            % shift clusters larger than deleted down by 1
            if numel(deleted) == 1
                gtMask = (spikeClusters > deleted);
                spikeClusters(gtMask) = spikeClusters(gtMask) - 1;
            else
                keepMe = setdiff(1:nClustersOld, deleted);
                nClustersNew = numel(keepMe);

                good = (spikeClusters > 0);
                mapFrom = zeros(1, nClustersOld);
                mapFrom(keepMe) = 1:nClustersNew;

                spikeClusters(good) = mapFrom(spikeClusters(good));
            end

        case 'merge'
            iCluster = op{3};
            jCluster = op{4};
            spikeClusters(spikeClusters == jCluster) = iCluster;

            % shift clusters larger than jCluster down by 1
            gtMask = (spikeClusters > jCluster);
            spikeClusters(gtMask) = spikeClusters(gtMask) - 1;

        case 'split'
            iCluster = op{3};
            retained = op{4};

            % shift clusters larger than iCluster up by 1 (make
            % room for splitted off cluster)
            gtMask = (spikeClusters > iCluster);
            spikeClusters(gtMask) = spikeClusters(gtMask) + 1;

            % take splitted off spikes and make a new cluster
            % of them
            iMask = (spikeClusters == iCluster);
            iSpikes = find(iMask);
            jSpikes = iSpikes(~ismember(iSpikes, retained)); % spikes to split off
            spikeClusters(jSpikes) = iCluster + 1;

        case 'partition'
            iCluster = op{3};
            assignPart = op{4};
            initialSpikes = find(spikeClusters == iCluster);

            for kSplit = 1:numel(assignPart)-1
                if kSplit > 1
                    iSpikes = find(spikeClusters == iCluster);
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
                mask = (spikeClusters > iCluster + kSplit - 1);
                spikeClusters(mask) = spikeClusters(mask) + 1;
                spikeClusters(splitOff) = iCluster + kSplit;
            end

        otherwise
            diffs = oldHistory{j, 3};
            if ~isempty(diffs)
                iDiffs = diffs(1, :);
                sDiffs = diffs(2, :);
                spikeClusters(iDiffs) = sDiffs;
            end

    end

    fwrite(fidHist, int32(j), 'int32');
    fwrite(fidHist, spikeClusters, 'int32');
    history(j) = op{2};
end

fclose(fidHist);
end