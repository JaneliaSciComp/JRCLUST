%--------------------------------------------------------------------------
function S_clu = S_clu_sort_(S_clu, criterion)
    % sort clusters by the centroid position
    % criterion: {'', 'clusterYPositions + clusterXPositions'}

    if nargin < 2
        criterion = '';
    end

    % Sort clusters by sites
    if isempty(criterion)
        criterion = 'clusterSites';
    end

    switch criterion
        case 'clusterYPositions + clusterXPositions'
            [~, viCluSort] = sort(S_clu.clusterYPositions + S_clu.clusterXPositions, 'ascend');
        otherwise
            [~, viCluSort] = sort(S_clu.(criterion), 'ascend');
    end

    S_clu.spikeClusters = mapIndex_(S_clu.spikeClusters, viCluSort);
    S_clu = struct_reorder_(S_clu, viCluSort, 'spikesByCluster', 'clusterXPositions', ...
        'clusterYPositions', 'vnSpk_clu', 'clusterSites', 'cviTime_clu', 'clusterNotes');

    S_clu = S_clu_refresh_(S_clu);
end %func
