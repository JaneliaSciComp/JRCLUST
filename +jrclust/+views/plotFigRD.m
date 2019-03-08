function hFigRD = plotFigRD(hFigRD, hClust, hCfg)
    %DOPLOTFIGRD Display rho-delta data   
    hFigRD.clf();

%     if isfield(hClust, 'cS_clu_shank')
%         cellfun(@(hc) doPlotFigRD(hc, P), hClust.cS_clu_shank);
%         return;
%     end

    if strcmp(hCfg.RDDetrendMode, 'global')
        [centers, x, y] = jrclust.sort.detrendRhoDelta(hClust, hClust.spikesBySite, 0, hCfg);
        y = jrclust.utils.nanlog10(y);
        fDetrend = 1;
    elseif strcmp(hCfg.RDDetrendMode, 'local')
        [centers, x, y] = jrclust.sort.detrendRhoDelta(hClust, hClust.spikesBySite, 1, hCfg);
        y = jrclust.utils.nanlog10(y);
        fDetrend = 1;
    else
        centers = find(hClust.spikeRho(:) > 10^(hCfg.log10RhoCut) & hClust.spikeDelta(:) > 10^(hCfg.log10DeltaCut));
        x = jrclust.utils.nanlog10(hClust.spikeRho(:));
        y = jrclust.utils.nanlog10(hClust.spikeDelta(:));
        fDetrend = 0;
    end

    indices = randsample(numel(x), min(numel(x), 50000), 0);

    % label cluster centers
    if ~isempty(hClust.clusterCenters)
        centers = hClust.clusterCenters; % do not overwrite
    end
    centersX = double(x(centers));
    centersY = double(y(centers));

    hFigRD.addPlot('allSpikes', x(indices), y(indices), '.');
    hFigRD.addPlot('centers', centersX, centersY, 'r.');
    hFigRD.axApply('default', @hold, 'on');

    hFigRD.axApply('default', @axis, 'tight');
    hFigRD.axApply('default', @axis, [-4 -.5 -1 2]);
    hFigRD.axApply('default', @set, 'XScale', 'linear', 'YScale', 'linear');

    % show rho/delta cutoff lines
    hFigRD.addPlot('RDCuts', hCfg.log10RhoCut * [1 1], hFigRD.axApply('default', @get, 'YLim'), 'r--', ...
                   hFigRD.axApply('default', @get, 'XLim'), hCfg.log10DeltaCut*[1, 1], 'r--');
    hFigRD.axApply('default', @grid, 'on');

    % set labels
    hFigRD.axApply('default', @xlabel, 'log10 rho');
    if fDetrend
        hFigRD.axApply('default', @ylabel, 'log10 delta (detrended)');
    else
        hFigRD.axApply('default', @ylabel, 'log10 delta');
    end
    hFigRD.axApply('default', @title, sprintf('rho-cut: %f, delta-cut: %f', hCfg.log10RhoCut, hCfg.log10DeltaCut));

    drawnow;
end