%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function flag = exist_file_(vcFile)
    % Different from exist(vcFile, 'file') which uses search path
    if isempty(vcFile)
        flag = 0;
    else
        flag = ~isempty(dir(vcFile));
    end
end %func
