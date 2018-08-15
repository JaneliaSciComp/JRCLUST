function twelve(filename)
    %TWELVE convert JRCv3-style prm files to JRCv4

    % WORK IN PROGRESS
    % TODO: use file2cellstr_, &c. (see updateParamFile.m)

    if ~endsWith(filename, '.prm')
        error('filename must be a .prm file');
    end

    % get the mapping between old and new names
    replaceMe = cell(0, 2);

    base = fileparts(fileparts(fullfile(mfilename('fullpath'))));
    fh = fopen(fullfile(base, 'params', 'prmRenames.txt'));

    tline = fgetl(fh);
    while ischar(tline)
        if ~isempty(tline)
            tline = strsplit(tline, '=');
            oldPrm = strtrim(tline{1});
            newPrm = strtrim(tline{2});

            if startsWith(oldPrm, 'P.')
                replaceMe(end+1, :) = {oldPrm(3:end), newPrm(3:end)};
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
            tline = strsplit(tline, '=');

            if numel(tline) > 2
                rhs = tline{2};
                for i=3:numel(tline)
                    rhs = [rhs '=' tline{i}];
                end
                tline = {tline{1}, rhs};
            end

            if numel(tline) == 2
                prm = strtrim(tline{1});
                prmVal = strtrim(tline{2});

                for i=1:size(replaceMe, 1)
                    replacePair = replaceMe(i, :);

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
        oldFilename = strrep(filename, '.prm', sprintf('-%s.prm', datestr(now,'yyyymmddTHHMMSS')));
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
    spkraw = strrep(filename, '.prm', '_spkraw.jrc');
    tracesBin = strrep(spkraw, 'spkraw.jrc', 'traces.bin');
    if isfile(spkraw)
        movefile(spkraw, tracesBin);
        fprintf('%s renamed to %s\n', spkraw, tracesBin);
    end

    spkwav = strrep(filename, '.prm', '_spkwav.jrc');
    waveformsBin = strrep(spkwav, 'spkwav.jrc', 'waveforms.bin');
    if isfile(spkwav)
        movefile(spkwav, waveformsBin);
        fprintf('%s renamed to %s\n', spkwav, waveformsBin);
    end

    spkfet = strrep(filename, '.prm', '_spkfet.jrc');
    featuresBin = strrep(spkfet, 'spkfet.jrc', 'features.bin');
    if isfile(spkfet)
        movefile(spkfet, featuresBin);
        fprintf('%s renamed to %s\n', spkfet, featuresBin);
    end
end
