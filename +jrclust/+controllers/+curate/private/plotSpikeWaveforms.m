function hFigWav = plotSpikeWaveforms(hFigWav, hClust, hCfg, maxAmp)
    %PLOTSPIKEWAVEFORMS Plot individual waveforms in the main view
    [xData, yData, showSites] = deal(cell(hClust.nClusters, 1));

    nSpikesCluster = zeros(hClust.nClusters, 1);
    siteNeighbors = hCfg.siteNeighbors(:, hClust.clusterSites);

    for iCluster = 1:hClust.nClusters
        try
            iSpikes = jrclust.utils.subsample(hClust.getCenteredSpikes(iCluster), hCfg.nSpk_show);
            iSites = siteNeighbors(:, iCluster);

            if hCfg.showRaw
                iWaveforms = jrclust.utils.rawTouV(hClust.spikesRaw(:, :, iSpikes), hCfg);
                iWaveforms = jrclust.filters.fftLowpass(iWaveforms, hCfg.getOr('fc_spkwav_show', []), hCfg.sampleRate);
            else
                iWaveforms = jrclust.utils.filtTouV(hClust.spikesFilt(:, :, iSpikes), hCfg);
            end

            [yData{iCluster}, xData{iCluster}] = wfToPlot(iWaveforms, iCluster, iSites, maxAmp, hCfg);
            showSites{iCluster} = iSites;
            nSpikesCluster(iCluster) = size(iWaveforms, 3);
        catch ME
            warning(ME.identifier, 'Can''t plot cluster %d: %s', iCluster, ME.message);
        end
    end

    userData = struct('showSites', showSites, 'nSpikesCluster', nSpikesCluster);
    if hFigWav.hasPlot('hSpkAll')
        hFigWav.updatePlot('hSpkAll', jrclust.utils.neCell2mat(xData), jrclust.utils.neCell2mat(yData), userData);
    else
        hFigWav.addPlot('hSpkAll', jrclust.utils.neCell2mat(xData), jrclust.utils.neCell2mat(yData), 'Color', [.5 .5 .5], 'LineWidth', .5, 'UserData', userData);
    end
end

%% LOCAL FUNCTIONS
function [yData, xData] = wfToPlot(waveforms, iCluster, iSites, maxAmp, hCfg)
    %WFTOPLOT Scale and translate waveforms by cluster ID and site
    iCluster = double(iCluster);

    if isempty(iSites)
        iSites = 1:size(waveforms, 2);
    end

    nSpikes = size(waveforms, 3);
    nSites = numel(iSites);

    waveforms = single(waveforms) / maxAmp;
    % orient waveforms in Y by which site they occur on
    waveforms = waveforms + repmat(single(iSites(:)'), [size(waveforms, 1), 1, size(waveforms, 3)]);
    waveforms([1, end], :, :) = nan;
    yData = waveforms(:);

    % x values are samples offset by cluster number
    xData = getXRange(iCluster, hCfg);
    xData = repmat(xData(:), [1, nSites * nSpikes]);
    xData = xData(:);
end

