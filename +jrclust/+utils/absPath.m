function ap = absPath(pathname, basedir)
    %ABSPATH Get the absolute path of pathname (file or directory)
    %   Returns empty str if pathname can't be found
    if nargin < 2 || ~exist(basedir, 'dir')
        basedir = pwd();
    elseif ~isAbsPath(basedir)
        basedir = fullfile(pwd(), basedir);
    end

    if nargin == 0 || ~ischar(pathname) || isempty(pathname)
        ap = '';
    elseif isAbsPath(pathname) && (exist(pathname, 'file') == 2 || exist(pathname, 'dir') == 7)
        ap = pathname;
    elseif any(pathname == '*') || any(pathname == '?')
        % first check full file
        d = dir(pathname);
        if ~isempty(d)
            dfolder = {d.folder};
            dname = {d.name};
            ap = arrayfun(@(i) fullfile(dfolder{i}, dname{i}), 1:numel(d), 'UniformOutput', 0); % cell array
        else % try with basedir hint
            d = dir(fullfile(basedir, pathname));
            ap = fullfile(basedir, {d.name});
            if isempty(ap)
                ap = '';
            end
        end
    else
        if exist(fullfile(basedir, pathname), 'file') == 2 || exist(fullfile(basedir, pathname), 'dir') == 7 % check hinted directory first
            ap = fullfile(basedir, pathname);
        elseif exist(fullfile(pwd(), pathname), 'file') || exist(fullfile(pwd(), pathname), 'dir') % try to find it relative to current directory
            ap = fullfile(pwd(), pathname);
        else % can't find it, you're on your own
            ap = '';
        end
    end
end

%% LOCAL FUNCTIONS
function iap = isAbsPath(pathname)
    %ISABSPATH Return true if pathname is an absolute path
    %   On *nix, it begins with /
    %   On Windows, it begins with / or \ after chopping off `[A-Za-z]:`
    if ispc() % Windows
        % chop off drive letter at beginning
        truncd = regexprep(pathname, '^[a-z]:', '', 'ignorecase');
        iap = ~strcmp(truncd, pathname) && (all(regexp(truncd, '^[/\\]')) == 1);
    else
        iap = regexp(pathname, '^/');
        iap = ~isempty(iap) && iap;
    end
end
