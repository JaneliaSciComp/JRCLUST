%--------------------------------------------------------------------------
% 10/8/17 JJJ: Created
function edit_(vcFile)
    % vcFile0 = vcFile;
    if ~exist_file_(vcFile)
        if matchFileExt_(vcFile, '.prb')
            vcFile = find_prb_(vcFile);
        else
            vcFile = jrcpath_(vcFile, 1);
        end
    end
    fprintf('Editing %s\n', vcFile);
    edit(vcFile);
end %func
