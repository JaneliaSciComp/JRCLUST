%--------------------------------------------------------------------------
function cursor_FigWavCor_(S0)
    if nargin==0, S0 = get(0, 'UseData'); end
    P = S0.P; S_clu = S0.S_clu;

    [hFig, S_fig] = get_fig_cache_('FigWavCor');
    if isempty(S_fig)
        [hFig, S_fig] = plot_FigWavCor_(S0);
    end
    iClu1 = S0.iCluCopy;
    if isempty(S0.iCluPaste)
        iClu2 = S0.iCluCopy;
    else
        iClu2 = S0.iCluPaste;
    end

    cor12 = S_clu.mrWavCor(iClu1, iClu2);
    set(S_fig.hCursorV, 'XData', iClu1*[1,1], 'YData', [.5, S_clu.nClu+.5]);
    title_(S_fig.hAx, sprintf('Clu%d vs. Clu%d: %0.3f; %s', iClu1, iClu2, cor12, S_fig.vcTitle));
    if iClu1==iClu2, color_H = [0 0 0]; else color_H = [1 0 0]; end
    set(S_fig.hCursorH, 'YData', iClu2*[1,1], 'XData', [.5, S_clu.nClu+.5], 'Color', color_H);
    xlim_(S_fig.hAx, trim_lim_(iClu1 + [-6,6], [.5, S_clu.nClu+.5]));
    ylim_(S_fig.hAx, trim_lim_(iClu2 + [-6,6], [.5, S_clu.nClu+.5]));
end %func
