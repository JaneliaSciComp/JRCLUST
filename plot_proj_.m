%--------------------------------------------------------------------------
function plot_proj_(hPlot, mrMin, mrMax, P, maxAmp)
    if nargin < 5
        [hFig, S_fig] = getCachedFig('FigProj');
        maxAmp = S_fig.maxAmp;
    end

    [vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, maxAmp, P.maxSite_show);

    % make struct
    maxPair = P.maxSite_show;
    sitesOfInterest = P.sitesOfInterest;
    S_plot = makeStruct_(mrMax, mrMin, sitesOfInterest, viPlot, tr_dim, maxPair, maxAmp);

    update_plot_(hPlot, vrX, vrY, S_plot);
end %func
