function twelve(filename)
    %TWELVE convert JRCv3-style prm files to JRCv4

    % WORK IN PROGRESS
    % TODO: use file2cellstr_, &c. (see updateParamFile.m)

    if ~endsWith(filename, '.prm')
        error('filename must be a .prm file');
    elseif ~exist(filename, 'file')
        error('%s not found', filename);
    end

    sessionName = strrep(filename, '.prm', '');
    timestamp = datestr(now,'yyyymmddTHHMMSS');

    % get the mapping between old and new names
    replaceMeP = cell(0, 2);
    replaceMeS0 = cell(0, 2);
    replaceMeS_clu = cell(0, 2);

    base = fileparts(fileparts(fullfile(mfilename('fullpath'))));
    fh = fopen(fullfile(base, 'params', 'prmRenames.txt'));

    tline = fgetl(fh);
    while ischar(tline)
        if ~isempty(tline)
            tline = strsplit(tline, '=');
            oldPrm = strtrim(tline{1});
            newPrm = strtrim(tline{2});

            if startsWith(oldPrm, 'P.')
                replaceMeP(end+1, :) = {oldPrm(3:end), newPrm(3:end)};
            elseif startsWith(oldPrm, 'S0')
                replaceMeS0(end+1, :) = {oldPrm(4:end), newPrm(4:end)};
            elseif startsWith(oldPrm, 'S_clu')
                replaceMeS_clu(end+1, :) = {oldPrm(7:end), newPrm(7:end)};
            end
        end

        tline = fgetl(fh);
    end
    fclose(fh);
    
    % write the new parameter file
    newlines = cell(0);
    fh = fopen(filename);
    nChanges = 0;

    tline = fgetl(fh);
    while ischar(tline)
        if isempty(tline)
            newlines{end+1} = '';
        else
            tline = strsplit(tline, ' = ');

            if numel(tline) > 2
                rhs = tline{2};
                for i=3:numel(tline)
                    rhs = [rhs ' = ' tline{i}];
                end
                tline = {tline{1}, rhs};
            end

            if numel(tline) == 2
                prm = strtrim(tline{1});
                prmVal = strtrim(tline{2});

                for i=1:size(replaceMeP, 1)
                    replacePair = replaceMeP(i, :);

                    if strcmp(prm, replacePair{1})
                        nChanges = nChanges + 1;
                        prm = replacePair{2};
                        break;
                    end
                end

                newlines{end+1} = sprintf('%s = %s', prm, prmVal);
            else
                newlines{end+1} = sprintf('%s', tline{1});
            end
        end
        tline = fgetl(fh);
    end
    fclose(fh);

    if nChanges > 0
        oldFilename = sprintf('%s-%s.prm', sessionName, timestamp);
        copyfile(filename, oldFilename);
        fprintf('Old parameter file has been saved to %s.\n', oldFilename);

        fh = fopen(filename, 'w');
        for iLine = 1:numel(newlines)
            fprintf(fh, '%s\n', newlines{iLine});
        end
        fclose(fh);
        fprintf('New parameter file written to %s.\n', filename);
    end

    % rename _spkraw.jrc, _spkwav.jrc, _spkfet.jrc
    spkraw = [sessionName '_spkraw.jrc'];
    tracesBin = strrep(spkraw, 'spkraw.jrc', 'traces.bin');
    if exist(spkraw, 'file')
        movefile(spkraw, tracesBin);
        fprintf('%s renamed to %s\n', spkraw, tracesBin);
    end

    spkwav = [sessionName '_spkwav.jrc'];
    waveformsBin = strrep(spkwav, 'spkwav.jrc', 'waveforms.bin');
    if exist(spkwav, 'file')
        movefile(spkwav, waveformsBin);
        fprintf('%s renamed to %s\n', spkwav, waveformsBin);
    end

    spkfet = [sessionName '_spkfet.jrc'];
    featuresBin = strrep(spkfet, 'spkfet.jrc', 'features.bin');
    if exist(spkfet, 'file')
        movefile(spkfet, featuresBin);
        fprintf('%s renamed to %s\n', spkfet, featuresBin);
    end

    % update variable names in the main .mat file
    matFile = [sessionName '_jrc.mat'];
    if exist(matFile, 'file')
        nChanges = 0;

        S0 = load(matFile, '-mat');
        % update S0
        for i = 1:size(replaceMeS0, 1)
            oldField = replaceMeS0{i, 1};
            newField = replaceMeS0{i, 2};
            if isfield(S0, oldField)
                nChanges = nChanges + 1;
                S0.(newField) = S0.(oldField);
                S0 = rmfield(S0, oldField);
            end
        end
        
        % update P
        P = S0.P;
        for i = 1:size(replaceMeP, 1)
            oldField = replaceMeP{i, 1};
            newField = replaceMeP{i, 2};
            if isfield(P, oldField)
                nChanges = nChanges + 1;
                P.(newField) = P.(oldField);
                P = rmfield(P, oldField);
            end
        end
        S0.P = P;
        
        % update S_clu
        S_clu = S0.S_clu;
        for i = 1:size(replaceMeS_clu, 1)
            oldField = replaceMeS_clu{i, 1};
            newField = replaceMeS_clu{i, 2};
            if isfield(S_clu, oldField)
                nChanges = nChanges + 1;
                S_clu.(newField) = S_clu.(oldField);
                S_clu = rmfield(S_clu, oldField);
            end
        end
        S0.S_clu = S_clu;
        
        if nChanges > 0
            % backup immediately
            oldFilename = sprintf('%s-%s_jrc.mat', sessionName, timestamp);
            copyfile(matFile, oldFilename);

            save(matFile, '-struct', 'S0');
            fprintf('Old MAT file has been saved to %s.\n', oldFilename);
        end
    end
    
    % update variable names in the _log.mat file
    logFile = [sessionName '_log.mat'];
    if exist(logFile, 'file')
        nChanges = 0;

        logStruct = load(logFile);
        % update log file
        for i = 1:size(replaceMeS_clu, 1)
            oldField = replaceMeS_clu{i, 1};
            newField = replaceMeS_clu{i, 2};
            if isfield(logStruct, oldField)
                nChanges = nChanges + 1;
                logStruct.(newField) = logStruct.(oldField);
                logStruct = rmfield(logStruct, oldField);
            end
        end

        if nChanges > 0
            % backup immediately
            oldFilename = sprintf('%s-%s_log.mat', sessionName, timestamp);
            copyfile(logFile, oldFilename);

            save(logFile, '-struct', 'logStruct');
            fprintf('Old log file has been saved to %s.\n', oldFilename);
        end
    end
end
