function metafile = findMeta(binfile)
    %FINDMETA Given a path to a recording (.bin, .dat), find its meta file
    %   Returns absolute path to meta file if it exists, otherwise empty
    %   string.
    [dirname, filename] = fileparts(binfile);
    metafile = jrclust.utils.absPath(fullfile(dirname, [filename '.meta']));
end

