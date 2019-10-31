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
        keepMe = int32([]);
        fidHist = fopen(obj.hCfg.histFile, 'r');

        % read in history to save...
        checkInt = fread(fidHist, 1, 'int32');
        while ~isempty(checkInt)
            if isKey(obj.history, checkInt)
                keepMe = [keepMe; checkInt; fread(fidHist, obj.nSpikes, 'int32')];
            else
                fseek(fidHist, 4*obj.nSpikes, 'cof');
            end
            checkInt = fread(fidHist, 1, 'int32');
        end
        fclose(fidHist);

        % ...and flush it back out to the hist file
        fidHist = fopen(obj.hCfg.histFile, 'w');
        fwrite(fidHist, keepMe, 'int32');
        fclose(fidHist);
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