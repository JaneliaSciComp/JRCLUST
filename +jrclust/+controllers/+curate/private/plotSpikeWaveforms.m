function hFigWav = plotSpikeWaveforms(hFigWav, hClust, hCfg, maxAmp)
    %PLOTSPIKEWAVEFORMS Plot individual waveforms in the main view
    [XData, YData, showSites] = deal(cell(hClust.nClusters, 1));

    nSpikesCluster = zeros(hClust.nClusters, 1);
    siteNeighbors = hCfg.siteNeighbors(:, hClust.clusterSites);

    for iCluster = 1:hClust.nClusters
        try
            iSpikes = jrclust.utils.subsample(hClust.getCenteredSpikes(iCluster), hCfg.nSpikesFigWav);
            iSites = siteNeighbors(:, iCluster);

            if hCfg.showRaw
                iWaveforms = jrclust.utils.rawTouV(hClust.spikesRaw(:, :, iSpikes), hCfg);
                iWaveforms = jrclust.filters.fftLowpass(iWaveforms, hCfg.getOr('fc_spkwav_show', []), hCfg.sampleRate);
            else
                iWaveforms = jrclust.utils.filtTouV(hClust.spikesFilt(:, :, iSpikes), hCfg);
            end

            [YData{iCluster}, XData{iCluster}] = wfToPlot(iWaveforms, iCluster, iSites, maxAmp, hCfg);
            showSites{iCluster} = iSites;
            nSpikesCluster(iCluster) = size(iWaveforms, 3);
        catch ME
            warning('Can''t plot cluster %d: %s', iCluster, ME.message);
        end
    end

    userData = struct('showSites', showSites, 'nSpikesCluster', nSpikesCluster);
    if hFigWav.hasPlot('hSpkAll')
        hFigWav.updatePlot('hSpkAll', jrclust.utils.neCell2mat(XData), jrclust.utils.neCell2mat(YData), userData);
    else
        hFigWav.addPlot('hSpkAll', jrclust.utils.neCell2mat(XData), jrclust.utils.neCell2mat(YData), 'Color', [.5 .5 .5], 'LineWidth', .5, 'UserData', userData);
    end
end

%% LOCAL FUNCTIONS
function [YData, XData] = wfToPlot(waveforms, iCluster, iSites, maxAmp, hCfg)
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
    YData = waveforms(:);

    % x values are samples offset by cluster number
    XData = getXRange(iCluster, hCfg);
    XData = repmat(XData(:), [1, nSites * nSpikes]);
    XData = XData(:);
end

