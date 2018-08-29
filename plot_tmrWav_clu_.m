%--------------------------------------------------------------------------
function [iClu, hPlot] = plot_tmrWav_clu_(S0, iClu, hPlot, vrColor)
    [S_clu, P] = getfield_(S0, 'S_clu', 'P');
    [hFig, S_fig] = getCachedFig('FigWav');
    if ~tryIsValid(hPlot)
        hPlot = plot(nan, nan, 'Color', vrColor, 'LineWidth', 2, 'Parent', S_fig.hAx);
    end
    if P.fWav_raw_show
        mrWav_clu1 = S_clu.tmrWav_raw_clu(:,:, iClu);
    else
        mrWav_clu1 = S_clu.tmrWav_clu(:,:, iClu);
    end
    multiplot(hPlot, S_fig.maxAmp, wav_clu_x_(iClu, P), mrWav_clu1);
    uistack_(hPlot, 'top');
end % function
