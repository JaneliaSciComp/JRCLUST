%--------------------------------------------------------------------------
function plotFeatureProjections(hPlot, mrMax, mrMin, P, maxAmp)
    if nargin < 5
        [~, figProjData] = getCachedFig('FigProj');
        maxAmp = figProjData.maxAmp;
    end

    switch lower(P.displayFeature)
        case {'vpp', 'vmin', 'vmax'}
            bounds = [0 maxAmp];
            maxPair = P.maxSite_show;
            vpp = 1;

        case 'kilosort'
            % round to nearest hundred on either side of 0
            bounds = round(max(abs([mrMin(:) ; mrMax(:)]))/100, 0)*[-100, 100];
            maxPair = [];
            vpp = 0;

        otherwise % not yet implemented
            error('display feature %s not yet implemented', P.displayFeature);

    end

    [xvals, yvals, viPlot, tr_dim] = featuresToSiteGrid(mrMax, mrMin, bounds, maxPair, vpp);

    % make struct
    maxPair = P.maxSite_show;
    sitesOfInterest = P.sitesOfInterest;
    S_plot = makeStruct_(mrMax, mrMin, sitesOfInterest, viPlot, tr_dim, maxPair, maxAmp);

    updatePlot(hPlot, xvals, yvals, S_plot);
end %func
