%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and testing
function csFiles = find_empty_files_(vcDir)
    % find files with 0 bytes

    if nargin==0, vcDir = []; end
    if isempty(vcDir), vcDir = pwd(); end
    vS_dir = dir(vcDir);
    viFile = find([vS_dir.bytes] == 0 & ~[vS_dir.isdir]);
    csFiles = {vS_dir(viFile).name};
    csFiles = cellfun(@(vc)[vcDir, filesep(), vc], csFiles, 'UniformOutput', 0);
end %func
