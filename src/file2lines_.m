%--------------------------------------------------------------------------
% Read a text file and output cell strings separated by new lines
% 7/24/17 JJJ: Code cleanup
function csLines = file2lines_(vcFile_file2struct)
    if ~exist(vcFile_file2struct, 'file')
        fprintf(2, '%s does not exist.\n', vcFile_file2struct);
        csLines = {};
        return;
    end

    fid = fopen(vcFile_file2struct, 'r');
    csLines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);

    csLines = csLines{1};
end %func
