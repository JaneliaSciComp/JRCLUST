function ap = absPath(filename, basedir)
    %ABSPATH Get the absolute path of filename
    %   Returns empty str if filename can't be found
    if nargin < 2
        basedir = pwd();
    end
    
    ap = '';

    % which returns empty if you're looking for a directory
    if isfile(filename) && isempty(which(filename))
        ap = filename;
    elseif isfile(filename)
        ap = which(filename);
    else % ~isfile(filename)
        filename = fullfile(basedir, filename);
        if isfile(filename) && isempty(which(filename)) % basedir is absolute
            ap = filename;
        elseif isfile(filename) % basedir is relative
            ap = which(filename);
        end % otherwise, can't find it, you're on your own
    end
end

