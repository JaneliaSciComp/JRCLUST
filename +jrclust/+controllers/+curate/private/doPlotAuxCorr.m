function doPlotAuxCorr(hClust, firingRates, auxSamples, auxChanCorr, auxTimes, selected)
    %DOPLOTAUXCORR Plot aux channel correlation figure
    hCfg = hClust.hCfg;
    [~, argsort] = sort(auxChanCorr, 'descend');

    nClustersShow = min(hCfg.nClustersShowAux, numel(auxChanCorr));
    auxLabel = hCfg.getOr('auxLabel', 'aux');
    nSubsamplesAux = hCfg.getOr('nSubsamplesAux', 100);

    if ~isempty(selected)
        nClustersShow = 1;
        argsort = selected;
    end

    hFigAux = jrclust.views.Figure('FigAux', [.5 0 .5 1], hCfg.sessionName, 1, 1);
    hTabGroup = hFigAux.figApply(@uitabgroup);

    for iiCluster = 1:nClustersShow
        iCluster = argsort(iiCluster);
        hTab = uitab(hTabGroup, 'Title', sprintf('Cluster %d', iCluster), 'BackgroundColor', 'w');
        axes('Parent', hTab);
        subplot(2, 1, 1);

        hAx = plotyy(auxTimes, firingRates(:, iCluster), auxTimes, auxSamples);

        xlabel('Time (s)');
        ylabel(hAx(1),'Firing Rate (Hz)');
        ylabel(hAx(2), auxLabel);

        iSite = hClust.clusterSites(iCluster);
        iTitle = sprintf('Cluster %d (Site %d, Chan %d): Corr=%0.3f', iCluster, iSite, hCfg.siteMap(iSite), auxChanCorr(iCluster));
        title(iTitle);
        set(hAx, 'XLim', auxTimes([1,end]));
        grid on;

        subplot(2, 1, 2);
        plot(auxSamples(1:nSubsamplesAux:end), firingRates(1:nSubsamplesAux:end,iCluster), 'k.');
        xlabel(auxLabel);
        ylabel('Firing Rate (Hz)');
        grid on;
    end
end
