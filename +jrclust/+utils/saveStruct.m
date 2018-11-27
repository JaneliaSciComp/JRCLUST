function saveStruct(S, filename)
    %SAVESTRUCT Save a struct to file in a robust way
    nRetry = 3;
    verYear = version('-release');
    verYear = str2double(verYear(1:end-1));

    if verYear >= 2017
        for iRetry = 1:nRetry
            try
                save(filename, '-struct', 'S', '-v7.3', '-nocompression');
                break;
            catch
                pause(.5);
            end

            warning('Could not save: %s', filename);
        end
    else
        for iRetry = 1:nRetry
            try
                save(filename, '-struct', 'S', '-v7.3');
                break;
            catch
                pause(.5);
            end

            warning('Could not save: %s', filename);
        end
    end
end
