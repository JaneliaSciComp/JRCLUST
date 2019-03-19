function commit(obj, msg)
    %COMMIT Commit a modification of clustering to history log
    %   Message format: parts are separated by a semicolon, no spaces
    %   1: datestr(timestamp, 31)
    %   2: operation (initial commit/import; delete; merge; split)
    %   3: cluster(s) operated on:  for delete, all clusters deleted, comma-separated
    %                               for merge, lower-valued cluster
    %                               for split, cluster splitted
    %   4: operation-specific data: for delete, empty
    %                               for merge, higher-valued cluster
    %                               for split, spikes retained, comma-separated
    %   Example: 2018-12-26 13:02:17;merge;4;5 denotes cluster 5
    %   merged into cluster 4
    if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
        warning('cannot branch from history; use revert() first');
        return;
    end

    % check for consistency before committing
    ic = obj.inconsistentFields();
    if ~isempty(ic)
        wmsg = strjoin(ic, '\n\t');
        warning('Cluster data inconsistent after previous operations:\n\t%s', wmsg);
        obj.editSeek(obj.nEdits);
        return;
    end

    % split commit message by commas
    msgParts = strsplit(lower(msg), ';');
    if numel(msgParts) < 2
        warning('malformed commit message: %s', msg);
        return
    end

    obj.history(end+1, 1:2) = msgParts(1:2);
    if startsWith(msgParts{2}, 'delete')
        % multiple clusters can be deleted at once
        deleted = cellfun(@(x) str2double(x), strsplit(msgParts{3}, ','));
        obj.history{end, 3} = deleted;
    elseif startsWith(msgParts{2}, 'merge')
        iCluster = str2double(msgParts{3});
        jCluster = str2double(msgParts{4});
        obj.history{end, 3} = min(iCluster, jCluster);
        obj.history{end, 4} = max(iCluster, jCluster);
    elseif startsWith(msgParts{2}, 'split')
        obj.history{end, 3} = str2double(msgParts{3});
        retained = cellfun(@(x) str2double(x), strsplit(msgParts{4}, ','));
        obj.history{end, 4} = retained;
    elseif startsWith(msgParts{2}, 'partition')
        obj.history{end, 3} = str2double(msgParts{3});
        obj.history{end, 4} = eval(strjoin(msgParts(4:end), ';')); % commit mechanism needs to change
    elseif any(obj.spikeClusters ~= obj.initialClustering) % store diff from previous clustering
        spikeClusters_ = obj.initialClustering;

        for j = 1:size(obj.history, 1)-1
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
        end % now caught up to the current state, can get diff

        iDiffs = find(obj.spikeClusters ~= spikeClusters_);
        sDiffs = obj.spikeClusters(iDiffs);
        obj.history{end, 3} = [iDiffs'; sDiffs'];
    end

    obj.editPos = size(obj.history, 1) - 1; % 0 for initial commmit, etc.
end