function [newFilename] = twelve(filename)
%TWELVE convert JRCv3-style prm files to JRCv4

    % WORK IN PROGRESS

    if ~endsWith(filename, '.prm')
        error('filename must be a .prm file');
    end
    
    newFilename = strrep(filename, '.prm', '-new.prm');
    
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
    oldFh = fopen(filename);
    newFh = fopen(newFilename, 'w');
    
    tline = fgetl(oldFh);
    while ischar(tline)
        if isempty(tline)
            fprintf(newFh, '\n');
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
                        prm = replacePair{2};
                        break;
                    end
                end

                fprintf(newFh, '%s = %s\n', prm, prmVal);
            else
                fprintf(newFh, '%s\n', tline{1});
            end
        end
        tline = fgetl(oldFh);
    end
    
    fclose(oldFh);
    fclose(newFh);
    
    fprintf(1, 'New parameter file written at %s\n', newFilename);
end

