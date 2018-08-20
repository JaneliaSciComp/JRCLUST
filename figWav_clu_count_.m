%--------------------------------------------------------------------------
function S_fig = figWav_clu_count_(S_fig, S_clu, fText)
    if nargin==0, [hFig, S_fig] = getCachedFig('FigWav'); end
    if nargin<2, S_clu = []; end
    if nargin<3, fText = 1; end
    if isempty(S_clu), S_clu = get0_('S_clu'); end

    if fText
        csText_clu = arrayfun(@(i)sprintf('%d(%d)', i, S_clu.vnSpk_clu(i)), 1:S_clu.nClusters, 'UniformOutput', 0);
    else
        csText_clu = arrayfun(@(i)sprintf('%d', i), 1:S_clu.nClusters, 'UniformOutput', 0);
    end
    set(S_fig.hAx, 'Xtick', 1:S_clu.nClusters, 'XTickLabel', csText_clu, 'FontSize', 8);
    try
        if fText
            xtickangle(S_fig.hAx, -20);
        else
            xtickangle(S_fig.hAx, 0);
        end
    catch;
    end

    S_fig.fText = fText;
    if nargout==0
        hFig = getCachedFig('FigWav');
        set(hFig, 'UserData', S_fig);
    end
end %func
