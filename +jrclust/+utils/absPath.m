function ap = absPath(pathname, basedir)
    %ABSPATH Get the absolute path of pathname (file or directory)
    % Returns empty str if pathname can't be found
    if nargin < 2 || ~isdir(basedir)
        basedir = pwd();
    elseif ~isAbsPath(basedir)
        basedir = fullfile(pwd(), basedir);
    end

    if isAbsPath(pathname) && (isfile(pathname) || isdir(pathname))
        ap = pathname;
    else
        if isfile(fullfile(basedir, pathname)) || isdir(fullfile(basedir, pathname)) % check hinted directory first
            ap = fullfile(basedir, pathname);
        elseif isfile(pathname) || isdir(pathname) % try to find it relative to current directory
            ap = fullfile(pwd(), pathname);
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