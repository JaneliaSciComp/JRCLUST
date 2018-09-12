%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function git_pull_(vcVersion)
    % https://github.com/drbenvincent/github-sync-matlab
    % startDir = cd();
    if nargin<1, vcVersion = ''; end

    repoURL = 'https://github.com/jamesjun/JRCLUST';
    % repoName = 'JRCLUST';
    try
        if isempty(vcVersion)
            code = system('git pull');
        else
            code = system(sprintf('git reset --hard "%s"', vcVersion));
        end
    catch
        code = -1;
    end
    if code ~= 0
        fprintf(2, 'Not a git repository. Please run the following command to clone from GitHub.\n');
        fprintf(2, '\tsystem(''git clone %s.git myDest''\n', repoURL);
        fprintf(2, '\tReplace "myDest" with the desired installation location or omit to install in ./JRCLUST.\n', repoURL);
        fprintf(2, '\tYou may need to install git from https://git-scm.com/downloads.\n');
    else
        edit('change_log.txt');
    end
end %func
