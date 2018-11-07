function verstr = version()
    %VERSION Get version info
    md = jrclust.utils.info;

    verstr = sprintf('%s %d.%d.%d', md.program, md.version.major, ...
                     md.version.minor, md.version.release);
end

