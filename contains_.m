%--------------------------------------------------------------------------
function flag = contains_(vc, cs)
    % check if vs contains any of cs
    if ischar(cs), cs = {cs}; end
    flag = 0;
    for i = 1:numel(cs)
        if ~isempty(strfind(vc, cs{i}))
            flag = 1;
            break;
        end
    end
end %func
