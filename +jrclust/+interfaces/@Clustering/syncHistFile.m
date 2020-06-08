function syncHistFile(obj)
    %SYNCHISTFILE Sync up history file with current history
    d = dir(obj.hCfg.histFile);

    if isempty(d) || obj.nSpikes == 0 % file does not exist
        fclose(fopen(obj.hCfg.histFile, 'w')); % truncate history file
        d = dir(obj.hCfg.histFile);
    end

    if obj.nSpikes == 0
        return;
    end

    nEntries = d.bytes / 4 / (obj.nSpikes + 1); % int32 spike table plus header

    if nEntries ~= ceil(nEntries) % non-integer number of "entries", likely corrupt file
        obj.history = resetHistFile(obj.hCfg.histFile, obj.initialClustering, obj.spikeClusters, obj.nEdits);
    elseif nEntries > obj.nEdits % more edits in file than in stored history, pare file down to match history
        keepMe = cell2mat(keys(obj.history));
        try
            mm = memmapfile(obj.hCfg.histFile, 'Format', {'int32', [obj.nSpikes + 1 nEntries], 'Data'}, 'Writable', true);
            editKeys = mm.Data.Data(1, :);
        catch ME
            if ispc
               %% hack for getting around windows path length limitation
               tmp=['\\?\',obj.hCfg.histFile];
               mm = memmapfile(tmp, 'Format', {'int32', [obj.nSpikes + 1 nEntries], 'Data'}, 'Writable', true);
               editKeys = mm.Data.Data(1, :);
            else
               rethrow(ME);
            end
        end
        iKeepMe = ismember(editKeys, keepMe);

        if ~jrclust.utils.isEqual(find(iKeepMe), 1:sum(iKeepMe))
            % take the subset of entries we want to keep...
            mm.Data.Data(:, 1:obj.nEdits) = mm.Data.Data(:, iKeepMe);
        end

        clear mm; % flush changes and close

        % ...and truncate the file after this
        java.io.RandomAccessFile(obj.hCfg.histFile, 'rw').setLength(4 * obj.nEdits * (obj.nSpikes + 1));
    elseif nEntries < obj.nEdits % more edits in stored history than in file, consolidate edits
        obj.history = resetHistFile(obj.hCfg.histFile, obj.initialClustering, obj.spikeClusters, obj.nEdits);
    elseif nEntries > 0
        fidHist = fopen(obj.hCfg.histFile, 'r');
        fRes = 0;
        while fRes > -1
            checkInt = fread(fidHist, 1, 'int32');
            if isempty(checkInt) % EOF
                fRes = -1;
                break
            elseif ~isKey(obj.history, checkInt)
                break;
            end
            fRes = fseek(fidHist, 4*obj.nSpikes, 'cof');
        end
        fclose(fidHist);

        if fRes == 0 % checkInt failed
            obj.history = resetHistFile(obj.hCfg.histFile, obj.initialClustering, obj.spikeClusters, obj.nEdits);
        end
    end

    obj.editPos = obj.nEdits;
end

%% LOCAL FUNCTIONS
function history = resetHistFile(histFile, initialClustering, spikeClusters, nEdits)
    history = containers.Map('KeyType', 'int32', 'ValueType', 'char');
    fidHist = fopen(histFile, 'w');

    if nEdits > 0
        % write initial clustering
        history(1) = 'initial commit';
        fwrite(fidHist, int32(1), 'int32');
        fwrite(fidHist, int32(initialClustering), 'int32');

        if any(spikeClusters ~= initialClustering)
            history(2) = 'update clustering';
            fwrite(fidHist, int32(2), 'int32');
            fwrite(fidHist, int32(spikeClusters), 'int32');
        end
    end % otherwise truncate file

    fclose(fidHist);
end
