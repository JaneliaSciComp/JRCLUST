function flag = is_legal_filename_(filename)
    %IS_LEGAL_FILENAME_ check whether a filename is licit

    % permitted characters are ASCII letters, numbers, hyphens,
    % underscores, and periods; filenames can't end with a period
    try
        [~, filename, ext] = fileparts(filename);
        flag = regexpi(filename, '^[a-z0-9\-_\.]+$') && ...
            (isempty(ext) || regexpi(ext, '^\.[a-z0-9\-_]+$'));
    catch
        flag = 0;
    end
end
