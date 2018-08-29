%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function csLines = file2cellstr_(vcFile)
    % read text file to a cell string

    fid = fopen(vcFile, 'r');
    csLines = {};
    while ~feof(fid), csLines{end+1} = fgetl(fid); end
    fclose(fid);
    csLines = csLines';
end % function
