function ap = absPath(filename, basedir)
    %ABSPATH Get the absolute path of filename
    % Returns empty str if filename can't be found
    if nargin < 2 || ~isdir(basedir)
        basedir = pwd();
    elseif ~isAbsPath(basedir)
        basedir = fullfile(pwd(), basedir);
    end

    if isAbsPath(filename) && isfile(filename)
        ap = filename;
    else
        if isfile(fullfile(basedir, filename)) % check hinted directory first
            ap = fullfile(basedir, filename);
        elseif isfile(filename) % try to find it relative to current directory
            ap = fullfile(pwd(), filename);
        else % can't find it, you're on your own
            ap = '';
        end
    end
end

%% LOCAL FUNCTIONS
function ia = isAbsPath(pathname)
    %ISABSPATH Return true if pathname is an absolute path
    % On *nix, it begins with /
    % On Windows, it begins with / or \ after chopping off `[A-Za-z]:`

    if ispc() % Windows
        % chop off drive letter at beginning
        truncd = regexprep(pathname, '^[a-z]:', '', 'ignorecase');
        ia = ~strcmp(truncd, pathname) && (all(regexp(truncd, '^[/\\]')) == 1);
    else
        ia = (all(regexp(pathname, '^/')) == 1);
    end
end