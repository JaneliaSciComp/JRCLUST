%--------------------------------------------------------------------------
function struct_save_(S, vcFile, fVerbose)
    % 7/13/17 JJJ: Version check routine
    nRetry = 3;
    if nargin<3, fVerbose = 0; end
    if fVerbose
        fprintf('Saving a struct to %s...\n', vcFile); t1=tic;
    end
    version_year = version('-release');
    version_year = str2double(version_year(1:end-1));
    if version_year >= 2017
        for iRetry=1:nRetry
            try
                save(vcFile, '-struct', 'S', '-v7.3', '-nocompression'); %faster
                break;
            catch
                pause(.5);
            end
            fprintf(2, 'Saving failed: %s\n', vcFile);
        end
    else
        for iRetry=1:nRetry
            try
                save(vcFile, '-struct', 'S', '-v7.3');
                break;
            catch
                pause(.5);
            end
            fprintf(2, 'Saving failed: %s\n', vcFile);
        end
    end
    if fVerbose
        fprintf('\ttook %0.1fs.\n', toc(t1));
    end
end %func
