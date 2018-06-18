%--------------------------------------------------------------------------
function plot_FigPos_(S_clu1, S_clu2)
    [hFig, S_fig] = get_fig_cache_('FigPos');
    [S0, P, S_clu] = get0_();

    % plot waveform in space
    if isempty(S_fig)
        S_fig.hAx = axes_new_(hFig);
    else
        cla(S_fig.hAx); hold(S_fig.hAx, 'on');
    end
    plot_unit_(S_clu1, S_fig.hAx, [0 0 0]);
    %vrPosXY1 = [S_clu.vrPosX_clu(S_clu1.iClu), S_clu.vrPosY_clu(S_clu1.iClu)] / P.um_per_pix;
    vrPosXY1 = [S_clu.vrPosX_clu(S_clu1.iClu), S_clu.vrPosY_clu(S_clu1.iClu)];
    nSpk1 = S_clu.vnSpk_clu(S_clu1.iClu);
    if isempty(S_clu2)
        vcTitle = sprintf('Unit %d: %d spikes; (X=%0.1f, Y=%0.1f) [um]', S_clu1.iClu, nSpk1, vrPosXY1);
        try
            vcTitle = sprintf('%s\n%0.1fuVmin, %0.1fuVpp, SNR:%0.1f IsoDist:%0.1f ISIrat:%0.2f L-rat:%0.1f', ...
            vcTitle, S_clu1.uVmin, S_clu1.uVpp, S_clu1.snr, S_clu1.iso_dist, S_clu1.isi_ratio, S_clu1.l_ratio);
        catch
        end
    else
        nSpk2 = S_clu.vnSpk_clu(S_clu2.iClu);
        vrPosXY2 = [S_clu.vrPosX_clu(S_clu2.iClu), S_clu.vrPosY_clu(S_clu2.iClu)] / P.um_per_pix;
        plot_unit_(S_clu2, S_fig.hAx, [1 0 0]);
        vcTitle = sprintf('Unit %d(black)/%d(red); (%d/%d) spikes\n(X=%0.1f/%0.1f, Y=%0.1f/%0.1f) [um]', ...
        S_clu1.iClu, S_clu2.iClu, nSpk1, nSpk2, ...
        [vrPosXY1(1), vrPosXY2(1), vrPosXY1(2), vrPosXY2(2)]);
    end
    title_(S_fig.hAx, vcTitle);
    set(hFig, 'UserData', S_fig);
end %func
