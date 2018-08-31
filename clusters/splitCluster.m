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
    % [S_clu, newCluster] = clu_reorder_(S_clu, oldCluster);
    newCluster = oldCluster + 1; % set to the cluster immediately after oldCluster
    if newCluster ~= S_clu.nClusters % old core of clu_reorder_
        % tag spikes affected by a move of the cluster at nClusters to position oldCluster + 1
        addOne = S_clu.spikeClusters > oldCluster & S_clu.spikeClusters < S_clu.nClusters;
        % all spikes assigned to last cluster now assigned to oldCluster + 1
        S_clu.spikeClusters(S_clu.spikeClusters == S_clu.nClusters) = newCluster;
        % increment cluster IDs of all tagged spikes
        S_clu.spikeClusters(addOne) = S_clu.spikeClusters(addOne) + 1;

        % remap data at index nClusters into index newCluster and shift the rest
        clusterRemap = [1:oldCluster, S_clu.nClusters, oldCluster+1:S_clu.nClusters-1];
        S_clu = S_clu_select_(S_clu, clusterRemap);
    end

    % update all the other views
    S_clu = updateSimScore(S_clu);
    [S_clu, S0] = S_clu_commit_(S_clu, 'splitCluster');
    plotFigWav(S0); % redraw plot
    plotFigClusterCor(S0);

    save_log_(sprintf('split %d', oldCluster), S0); %@TODO: specify which cut to use

    % select two clusters being split
    button_CluWav_simulate_(oldCluster, newCluster);
    tryClose(hMsg);

    fprintf('%s [W] splitted Clu %d\n', datestr(now, 'HH:MM:SS'), oldCluster);
    figure_wait_(0);
end % function
