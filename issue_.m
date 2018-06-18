%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created
function issue_(vcMode)
    switch lower(vcMode)
        case 'post', web('https://github.com/jamesjun/JRCLUST/issues/new')
        case 'search', web('https://github.com/jamesjun/JRCLUST/issues')
    end %switch
end %func
