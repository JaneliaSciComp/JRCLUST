function dirname = basedir()
    %BASEDIR Get root directory of JRCLUST instance
    dirname = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

