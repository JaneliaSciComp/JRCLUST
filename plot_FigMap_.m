%--------------------------------------------------------------------------
function plot_FigMap_(S0)
    if nargin<1, S0 = get(0, 'UserData'); end
    P = S0.P; S_clu = S0.S_clu;
    [hFig, S_fig] = getCachedFig('FigMap');

    mrWav1 = S_clu.tmrWav_clu(:,:,S0.iCluCopy);
    vrVpp = squeeze_(max(mrWav1) - min(mrWav1));
    mrVpp = repmat(vrVpp(:)', [4, 1]);
    if isempty(S_fig)
        S_fig.hAx = axes_new_(hFig);
        [S_fig.mrPatchX, S_fig.mrPatchY] = probe_map_(P);
        S_fig.hPatch = patch(S_fig.mrPatchX, S_fig.mrPatchY, mrVpp, ...
        'EdgeColor', 'k', 'FaceColor', 'flat');
        S_fig.alim = [min(S_fig.mrPatchX(:)), max(S_fig.mrPatchX(:)), min(S_fig.mrPatchY(:)), max(S_fig.mrPatchY(:))];
        colormap jet;
        mouse_figure(hFig);
        nSites = size(P.mrSiteXY,1);
        csText = arrayfun(@(i)sprintf('%d', i), 1:nSites, 'UniformOutput', 0);
        S_fig.hText = text(P.mrSiteXY(:,1), P.mrSiteXY(:,2), csText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        xlabel('X Position (\mum)');
        ylabel('Y Position (\mum)');
    else
        set(S_fig.hPatch, 'CData', mrVpp);
    end
    title_(S_fig.hAx, sprintf('Max: %0.1f \\muVpp', max(vrVpp)));
    axis_(S_fig.hAx, S_fig.alim);
    caxis(S_fig.hAx, [0, max(vrVpp)]);

    set(hFig, 'UserData', S_fig);
end %func
