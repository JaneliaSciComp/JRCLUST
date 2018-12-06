function hFigRD = doPlotFigRD(hFigRD, hClust, hCfg)
    %DOPLOTFIGRD Display rho-delta data   
    hFigRD.clf();

%     if isfield(hClust, 'cS_clu_shank')
%         cellfun(@(hc) doPlotFigRD(hc, P), hClust.cS_clu_shank);
%         return;
%     end

    if strcmp(hCfg.rlDetrendMode, 'none')
        centers = find(hClust.spikeRho(:) > 10^(hCfg.log10RhoCut) & hClust.spikeDelta(:) > 10^(hCfg.log10DeltaCut));
        x = jrclust.utils.nanlog10(hClust.spikeRho(:));
        y = jrclust.utils.nanlog10(hClust.spikeDelta(:));
        fDetrend = false;
    elseif strcmp(hCfg.rlDetrendMode, 'global')
        [centers, x, y] = jrclust.clustering.detrendRhoDelta(hClust, hClust.spikesBySite, false, hCfg);
        y = jrclust.utils.nanlog10(y);
        fDetrend = true;
    elseif strcmp(hCfg.rlDetrendMode, 'local')
        [centers, x, y] = jrclust.clustering.detrendRhoDelta(hClust, hClust.spikesBySite, true, hCfg);
        y = jrclust.utils.nanlog10(y);
        fDetrend = true;
    end

    hFigRD.addPlot('allSpikes', x, y, '.');
    hFigRD.hold('on');

    hFigRD.axis('tight');
    hFigRD.axis([-4 -.5 -1 2]);
    hFigRD.axSet('XScale', 'linear', 'YScale', 'linear');

    % show rho/delta cutoff lines
    hFigRD.addPlot('RDCuts', hCfg.log10RhoCut * [1 1], hFigRD.axGet('YLim'), 'r--', ...
                   hFigRD.axGet('XLim'), hCfg.log10DeltaCut*[1, 1], 'r--');
    hFigRD.grid('on');

    hFigRD.figData = struct('x', x, 'y', y);

    % label cluster centers
    if ~isempty(hClust.clusterCenters)
        centers = hClust.clusterCenters; % do not overwrite
    end
    centersX = double(x(centers));
    centersY = double(y(centers));
    hFigRD.addPlot('centers', centersX, centersY, 'r.');

    % set labels
    hFigRD.xlabel('log10 rho');
    if fDetrend
        hFigRD.ylabel('log10 delta (detrended)');
    else
        hFigRD.ylabel('log10 delta');
    end
    hFigRD.title(sprintf('rho-cut: %f, delta-cut: %f', hCfg.log10RhoCut, hCfg.log10DeltaCut));

    drawnow;
end
