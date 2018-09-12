%--------------------------------------------------------------------------
function nBytes = getBytes_(vcFile)
    S_dir = dir(vcFile);
    if isempty(S_dir), nBytes=[]; return; end
    nBytes = S_dir(1).bytes;
end %func
