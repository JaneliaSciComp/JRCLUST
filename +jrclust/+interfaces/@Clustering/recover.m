function flag = recover(obj, ask)
    %RECOVER Inconsistent data recovery, with a hammer
    %returns -1 if canceled by user,
    %         0 if failed
    %         1 if succeeded and spikeClusters DID NOT change
    %         2 if succeeded and spikeClusters DID change
    if nargin < 2
        ask = 0;
    end

    flag = 1;

    if isempty(obj.inconsistentFields())
        return;
    end

    spikeClusters = obj.spikeClusters;
    goodClusters = spikeClusters(spikeClusters > 0);

    nClusters = max(goodClusters);

    % first ensure all clusters are contiguous
    missingClusters = setdiff(1:nClusters, goodClusters);
    if ~isempty(missingClusters) && ask
        dlgAns = questdlg('Gaps found in spike table. Should I correct this?', 'Confirm fix');
        if ismember(dlgAns, {'No', 'Cancel', ''})
            flag = -1; % canceled by user
            return;
        end
    end

    for i = numel(missingClusters):-1:1
        iC = missingClusters(i);
        mask = spikeClusters > iC;
        spikeClusters(mask) = spikeClusters(mask) - 1;
    end

    nClusters = nClusters - numel(missingClusters);

    % confirm this is the expected number of clusters
    if ask
        question = sprintf('%d clusters found. Does this look correct?', nClusters);
        dlgAns = questdlg(question, 'Confirm unit count');

        if ismember(dlgAns, {'No', 'Cancel', ''})
            flag = -1; % canceled by user
            return;
        end
    end

    % backup existing fields
    backup = struct();
    if any(spikeClusters ~= obj.spikeClusters)
        backup.spikeClusters = obj.spikeClusters;
        obj.spikeClusters = spikeClusters;
    end

    inconsistentFields = obj.inconsistentFields();
    for i = 1:numel(inconsistentFields)
        pn = inconsistentFields{i};
        backup.(pn) = obj.(pn);
    end

    try
        if isfield(backup, 'spikeClusters') || ismember('spikesByCluster', inconsistentFields)
            % recompute spikesByCluster, unitFields, &c.
            obj.spikesByCluster = arrayfun(@(iC) find(spikeClusters == iC), 1:nClusters, 'UniformOutput', 0);
            
            % update count of spikes per unit
            obj.unitCount = cellfun(@numel, obj.spikesByCluster);

            % update cluster sites
            obj.clusterSites = double(cellfun(@(sC) mode(obj.spikeSites(sC)), obj.spikesByCluster));
        end

        % any of the waveform fields inconsistent?
        if isfield(backup, 'spikeClusters') || ~cellfun(@isempty, regexp(inconsistentFields, '.+Wf.+'))
            % update mean waveforms for all units
            obj.updateWaveforms([]);

            % update unit positions
            obj.computeCentroids([]);

            % compute quality scores for all units
            obj.computeQualityScores([]);
        end

        % correct notes if we need to
        if numel(obj.clusterNotes) < obj.nClusters
            obj.clusterNotes = [obj.clusterNotes(:); repmat({''}, obj.nClusters - numel(obj.clusterNotes), 1)];
        elseif numel(obj.clusterNotes) > obj.nClusters
            % don't remove any nonempty notes
            emptyCells = find(cellfun(@isempty, obj.clusterNotes));
            removeMe = emptyCells(end-numel(obj.clusterNotes)+obj.nClusters+1:end);
            
            % any nonempty notes getting shifted down?
            shiftingDown = find(removeMe < numel(obj.clusterNotes) - numel(removeMe));
            if ~isempty(shiftingDown)
                oldPositions = removeMe(shiftingDown) + 1;
                newPositions = oldPositions - [1; diff(oldPositions)];
                warning('notes at for units %s have now been shifted to %s. You should double check that these notes apply to the correct units.', ...
                    strjoin(arrayfun(@num2str, oldPositions, 'UniformOutput', 0), ','), strjoin(arrayfun(@num2str, newPositions, 'UniformOutput', 0), ','));
            end
            obj.clusterNotes(removeMe) = [];
        end
        
        flag = 1 + isfield(backup, 'spikeClusters');
    catch ME
        warning('restore failed, reverting');
        backedUpFields = fieldnames(backup);

        for i = 1:numel(backedUpFields)
            fn = backedUpFields{i};
            obj.(fn) = backup.(fn);
        end

        flag = 0;
    end
end

