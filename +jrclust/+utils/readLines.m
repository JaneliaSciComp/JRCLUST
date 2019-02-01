function lines = readLines(filename)
    %READLINES Read the lines of a text file into a cell array
    assert(exist(filename, 'file') == 2, '''%s'' is not a file', filename);

    fid = fopen(filename, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);

    lines = lines{1};
end

