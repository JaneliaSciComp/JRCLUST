%--------------------------------------------------------------------------
% 9/28/17 JJJ: Created and tested
function wiki_download_()
    repoURL = 'https://github.com/jamesjun/JRCLUST.wiki.git';
    repoName = 'wiki';
    if isempty(dir(repoName))
        vcEval = sprintf('git clone %s %s', repoURL, repoName);
    else
        vcEval = sprintf('git pull %s', repoName);
    end
    try
        code = system(vcEval);
    catch
        code = -1;
    end
    if code == 0
        fprintf('Wiki on GitHub is downloaded to ./%s/\n', repoName);
    else
        fprintf(2, 'Please install git from https://git-scm.com/downloads.\n');
    end
end %func
