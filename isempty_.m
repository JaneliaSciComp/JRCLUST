%--------------------------------------------------------------------------
function vl = isempty_(cvr)
    if iscell(cvr)
        vl = cellfun(@isempty, cvr);
    else
        vl = isempty(cvr);
    end
end % function
