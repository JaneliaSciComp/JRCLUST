%--------------------------------------------------------------------------
function S_clu = S_clu_sort_(S_clu, vcField_sort)
    % sort clusters by the centroid position
    % vcField_sort: {'', 'vrPosY_clu + vrPosX_clu'}

    if nargin<2, vcField_sort = ''; end

    % Sort clusters by its sites
    if isempty(vcField_sort), vcField_sort = 'clusterSites'; end

    switch vcField_sort
        case 'vrPosY_clu + vrPosX_clu'
        [~, viCluSort] = sort(S_clu.vrPosY_clu + S_clu.vrPosX_clu, 'ascend');
        otherwise
        [~, viCluSort] = sort(S_clu.(vcField_sort), 'ascend');
    end
    S_clu.spikeClusters = mapIndex_(S_clu.spikeClusters, viCluSort);
    S_clu = struct_reorder_(S_clu, viCluSort, ...
    'cviSpk_clu', 'vrPosX_clu', 'vrPosY_clu', 'vnSpk_clu', 'clusterSites', 'cviTime_clu', 'csNote_clu');
    % if isfield(S_clu, 'tmrWav_clu')
    %     S_clu.tmrWav_clu = S_clu.tmrWav_clu(:, :, viCluSort);
    % end
    S_clu = S_clu_refresh_(S_clu);
end %func
