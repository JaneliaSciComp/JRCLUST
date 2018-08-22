%--------------------------------------------------------------------------
function S_clu = splitCluster(oldCluster, spikesToSplitOff)
    % split cluster
    % update struct elements to reflect splitting

    figure_wait_(1);
    drawnow;

    [P, S_clu, spikeSites] = get0_('P', 'S_clu', 'spikeSites');
    hMsg = msgbox_open_('Splitting...');
    figure(getCachedFig('FigWav'));

    % create a new cluster (add at the end)
    nSpikesNew = sum(spikesToSplitOff); % number of spikes to split off
    newCluster = max(S_clu.spikeClusters) + 1;

    % update cluster count and index
    S_clu.nClusters = double(newCluster);
    S_clu.nSpikesPerCluster(oldCluster) = S_clu.nSpikesPerCluster(oldCluster) - nSpikesNew;
    S_clu.nSpikesPerCluster(newCluster) = sum(spikesToSplitOff);

    oldSpikes = find(S_clu.spikeClusters==oldCluster);
    newSpikes = oldSpikes(spikesToSplitOff);
    oldSpikes = oldSpikes(~spikesToSplitOff);

    % assign sites to old and new clusters
    oldSite = mode(spikeSites(oldSpikes));
    newSite = mode(spikeSites(newSpikes));

    if oldSite > newSite % swap order if old cluster is on a higher site
        [newSite, oldSite] = deal(oldSite, newSite);
        [newSpikes, oldSpikes] = deal(oldSpikes, newSpikes);
        spikesToSplitOff = ~spikesToSplitOff;
    end

    % update spikesByCluster and clusterSites
    S_clu.spikesByCluster{oldCluster} = oldSpikes;
    S_clu.clusterSites(oldCluster) = oldSite;

    S_clu.spikesByCluster{newCluster} = newSpikes;
    S_clu.clusterSites(newCluster) = newSite;

    % erase annotation on old cluster
    % try
    S_clu.clusterNotes{oldCluster} = '';
    S_clu.clusterNotes{end+1} = ''; %add another entry
    % catch
    % end

    % update spikeClusters for newSpikes
    S_clu.spikeClusters(newSpikes) = newCluster;
    S_clu = S_clu_update_(S_clu, [oldCluster, newCluster], P);

    % Bring the new cluster right next to the old one using index swap
    [S_clu, newCluster] = clu_reorder_(S_clu, oldCluster);

    % update all the other views
    [S_clu, S0] = S_clu_commit_(S_clu, 'splitCluster');
    plotFigWav(S0); % redraw plot
    plotFigClusterCor(S0);

    save_log_(sprintf('split %d', oldCluster), S0); %@TODO: specify which cut to use

    % select two clusters being split
    button_CluWav_simulate_(oldCluster, newCluster);
    tryClose(hMsg);

    fprintf('%s [W] splitted Clu %d\n', datestr(now, 'HH:MM:SS'), oldCluster);
    figure_wait_(0);
end %func
