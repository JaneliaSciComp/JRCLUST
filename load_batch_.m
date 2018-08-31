%--------------------------------------------------------------------------
% 8/4/17 JJJ: Function created, tested and documented
function csLines = load_batch_(vcFile_batch)
    % Load list of files in the text file, strip comments and empty lines
    % Remove empty lines, comment lines starting with '%' character

    csLines = file2cellstr_(vcFile_batch);
    % csFiles_prm = importdata(vcFile_batch);

    csLines = csLines(cellfun(@(x)~isempty(x), csLines)); % remove empty lines
    csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0); % remove empty space
    csLines = csLines(cellfun(@(x)~isempty(x), csLines)); %remove empty lines
    csLines = csLines(cellfun(@(x)x(1)~='%', csLines)); % remove comment lines

    % Removing comments that starts with "%"
    % func_comment = @(vc)vc(1) == '%';
    % viComment = cellfun(@(vc)func_comment(strtrim(vc)), csFiles_prm);
    % csFiles_prm(viComment) = [];
end % function
