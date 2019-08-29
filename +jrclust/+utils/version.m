function verstr = version()
    %VERSION Get version info
    md = jrclust.utils.info();
    verstr = sprintf('%d.%d.%d', md.version.major, md.version.minor, md.version.release);
    if ~isempty(md.version.tag)
        verstr = sprintf('%s-%s', verstr, md.version.tag);
    end

    if ~isempty(md.version.codename)
        verstr = sprintf('%s "%s"', verstr, md.version.codename);
    end

    if isfield(md.version, 'commitHash')
        idx = min(7, numel(md.version.commitHash));
        verstr = sprintf('%s (%s)', verstr, md.version.commitHash(1:idx));
    end
end

