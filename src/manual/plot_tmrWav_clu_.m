%--------------------------------------------------------------------------
function [iClu, hPlot] = plot_tmrWav_clu_(S0, iClu, hPlot, vrColor)
    [hClust, hCfg] = getfield_(S0, 'S_clu', 'P');
    [hFig, S_fig] = get_fig_cache_('FigWav');
    if ~isvalid_(hPlot)
        hPlot = plot(nan, nan, 'Color', vrColor, 'LineWidth', 2, 'Parent', S_fig.hAx);
    end
    if hCfg.fWav_raw_show
        mrWav_clu1 = hClust.tmrWav_raw_clu(:,:, iClu);
    else
        mrWav_clu1 = hClust.tmrWav_clu(:,:, iClu);
    end
    multiplot(hPlot, S_fig.maxAmp, getXRange(iClu, hCfg), mrWav_clu1);
    uistack_(hPlot, 'top');
end %func
