%--------------------------------------------------------------------------
function plotFeatureProjections(hPlot, mrMin, mrMax, P, maxAmp)
    if nargin < 5
        [hFig, S_fig] = getCachedFig('FigProj');
        maxAmp = S_fig.maxAmp;
    end

    switch lower(P.displayFeature)
        case 'kilosort'
            % round to nearest hundred on either side of 0
            bounds = round(max(abs([mrMin(:) ; mrMax(:)]))/100, 0)*[-100, 100];
            maxPair = [];

        otherwise % vpp et al.
            bounds = [0 maxAmp];
            maxPair = P.maxSite_show;
    end
    
    [vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, bounds, maxPair);

    % make struct
    maxPair = P.maxSite_show;
    sitesOfInterest = P.sitesOfInterest;
    S_plot = makeStruct_(mrMax, mrMin, sitesOfInterest, viPlot, tr_dim, maxPair, maxAmp);

    update_plot_(hPlot, vrX, vrY, S_plot);
end %func
