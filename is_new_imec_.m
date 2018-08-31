%--------------------------------------------------------------------------
function flag = is_new_imec_(vcFile)
    flag = ~isempty(regexp(lower(vcFile), '.imec.ap.bin$'));
end % function
