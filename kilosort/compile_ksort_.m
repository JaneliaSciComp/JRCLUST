%--------------------------------------------------------------------------
% Compile Kilosort code
% 10/5/17 JJJ: Error messages converted to warning
% 7/26/17 JJJ: Code cleanup and test
function fSuccess = compile_ksort_()
    oldwd = pwd();
    [dirname, ~] = fileparts(fullfile(mfilename('fullpath')));
    cd(fullfile(dirname, 'KiloSort', 'CUDA'));
    
    fSuccess = 1;
    nTry = 3;
    csFiles_cu = {'mexMPmuFEAT.cu', 'mexWtW2.cu', 'mexMPregMU.cu'};
    %delete ./kilosort/*.mex*;
    %delete ./kilosort/*.lib;
    for iFile = 1:numel(csFiles_cu)
        for iTry = 1:nTry
            try
                drawnow;
                eval(sprintf('mexcuda -largeArrayDims -v %s;', csFiles_cu{iFile}));
                fprintf('Kilosort compile success for %s.\n', csFiles_cu{iFile});
                break;
            catch
                if iTry == nTry
                    fprintf('\tKilosort could not be compiled: %s\n', csFiles_cu{iFile});
                    fSuccess = 0;
                end
            end
        end
    end
    if ~fSuccess
        fprintf('\tWarning: Kilosort could not be compiled but it may work fine. If not, install Visual Studio 2013 and run "jrc install".\n');
    end
    
    cd(oldwd);
end %func
