%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function flag = fileExists(filename)
    % Different from exist(filename, 'file') which uses search path
    if isempty(filename)
        flag = 0;
    else
        flag = ~isempty(dir(filename));
    end
end % function
