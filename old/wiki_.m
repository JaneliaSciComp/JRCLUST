%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created
function wiki_(vcPage)
    if nargin<1, vcPage = ''; end
    if isempty(vcPage)
        web('https://github.com/jamesjun/JRCLUST/wiki');
    else
        web(['https://github.com/jamesjun/JRCLUST/wiki/', vcPage]);
    end
end
