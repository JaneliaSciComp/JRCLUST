%--------------------------------------------------------------------------
function plotFeatureProjections(hPlot, mrMax, mrMin, P, maxAmp)
    if nargin < 5
        [~, figProjData] = getCachedFig('FigProj');
        maxAmp = figProjData.maxAmp;
    end

    switch lower(P.displayFeature)
        case {'vpp', 'vmin', 'vmax'}
            bounds = maxAmp*[0 1];
            maxPair = P.maxSite_show;
            vpp = 1;

            plotData = struct('mrMax', mrMax, 'mrMin', mrMin);

        case 'kilosort'
            % round up to nearest 50 on either side of 0
            bounds = maxAmp*[-1 1];
            maxPair = [];
            vpp = 0;

            plotData = struct('PC1', mrMax, 'PC2', mrMin);

        otherwise % not yet implemented
            error('display feature %s not yet implemented', P.displayFeature);

    end

    [xvals, yvals, viPlot, tr_dim] = featuresToSiteGrid(mrMax, mrMin, bounds, maxPair, vpp);

    % make struct
    maxPair = P.maxSite_show;
    sitesOfInterest = P.sitesOfInterest;

    plotData.sitesOfInterest = sitesOfInterest;
    plotData.viPlot = viPlot;
    plotData.tr_dim = tr_dim;
    plotData.maxPair = maxPair;
    plotData.maxAmp = maxAmp;
    % figData = makeStruct(mrMax, mrMin, sitesOfInterest, viPlot, tr_dim, maxPair, maxAmp);

    updatePlot(hPlot, xvals, yvals, plotData);
end % function
