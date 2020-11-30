function success = saveStruct(S, filename)
%SAVESTRUCT Save a struct to file in a somewhat more robust way.
success = false;

nRetries = 3;
verYear = version('-release');
verYear = str2double(verYear(1:end-1));

if verYear >= 2017
    for iRetry = 1:nRetries
        try
            save(filename, '-struct', 'S', '-v7.3', '-nocompression');
            success = true;
            break;
        catch ME
            pause(.5);
        end

        if ispc
            [success, msg] = saveTemp
            AndMove(S, filename, '-nocompression');
            if success
                break;
            end

            warning('JRC:saveStruct:noSave', 'Could not save %s: %s', filename, msg);
        else
            warning('JRC:saveStruct:noSave', 'Could not save %s', filename);                
        end
    end
else
    for iRetry = 1:nRetries
        try
            save(filename, '-struct', 'S', '-v7.3');
            success = true;
            break;
        catch
            pause(.5);
        end

        if ispc
            [success, msg] = saveTempAndMove(S, filename);
            if success
                break;
            end

            warning('JRC:saveStruct:noSave', 'Could not save %s: %s', filename, msg);
        else
            warning('JRC:saveStruct:noSave', 'Could not save %s', filename);
        end
    end
end
end %fun

%% LOCAL FUNCTION
function [success, msg] = saveTempAndMove(S, filename, varargin)
%SAVETEMPANDMOVE Try to save struct to a temporary location and move it to
%its intended destination.
%   Workaround for Windows path length limitations.
success = false;
tempResFile = [tempname() '.mat'];
save(tempResFile, '-struct', 'S', '-v7.3', varargin{:});
if exist(tempResFile, 'file') == 2
    [success, msg] = copyfile(tempResFile, ['\\?\' filename]);
end
end %fun
