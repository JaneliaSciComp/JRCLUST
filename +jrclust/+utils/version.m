function verstr = version()
    %VERSION Get version info
    md = jrclust.utils.info;
    verstr = sprintf('%d.%d.%d', md.version.major, md.version.minor, md.version.release);
end

