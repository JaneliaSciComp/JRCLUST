function saveStruct(S, filename) %#ok<INUSL>
    %SAVESTRUCT Save a struct to file in a robust way
    nRetries = 3;
    verYear = version('-release');
    verYear = str2double(verYear(1:end-1));

    if verYear >= 2017
        for iRetry = 1:nRetries
            try
                save(filename, '-struct', 'S', '-v7.3', '-nocompression');
                break;
            catch
                pause(.5);
            end
            if ispc 
               userDir = char(java.lang.System.getProperty('user.home'));
               filename_temp = [userDir,filesep,'tmp_jrclust.mat']; 
               save(filename_temp, '-struct', 'S', '-v7.3', '-nocompression');
               failed = system(sprintf('move "%s" "%s"',filename_temp,filename));
               if failed
                    warning('Could not save: %s', filename);
               end  
            else
                warning('Could not save: %s', filename);                
            end                
        end
    else
        for iRetry = 1:nRetries
            try
                save(filename, '-struct', 'S', '-v7.3');
                break;
            catch
                pause(.5);
            end
            if ispc 
               userDir = char(java.lang.System.getProperty('user.home'));
               filename_temp = [userDir,filesep,'tmp.mat']; 
               save(filename_temp, '-struct', 'S', '-v7.3');
               failed = system(sprintf('move "%s" "%s"',filename_temp,filename));
               if failed
                    warning('Could not save: %s', filename);
               end  
            else
                warning
        end
    end
end
