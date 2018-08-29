%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation. Added strtrim
function cs = first_string_(cs)
    % Return the first string, which is typically a variable name

    if ischar(cs), cs = {cs}; end

    for i=1:numel(cs)
        cs{i} = strtrim(cs{i});
        if isempty(cs{i}), continue; end
        cs1 = textscan(cs{i}, '%s', 'Delimiter', {' ','='});
        cs1 = cs1{1};
        cs{i} = cs1{1};
    end
end % function
