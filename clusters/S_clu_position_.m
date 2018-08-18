%--------------------------------------------------------------------------
function S_clu = S_clu_position_(S_clu, viClu_update)
    % determine cluster position from spike position
    % 6/27/17 JJJ: multiple features supported (single dimension such as energy and Vpp)
    global spikeFeatures

    if isempty(spikeFeatures)
        spikeFeatures = get_spkfet_(P);
    end

    if nargin<2, viClu_update = []; end
    P = get0_('P'); %P = S_clu.P;
    if ~isfield(S_clu, 'clusterXPositions'), S_clu.clusterXPositions = []; end
    if ~isfield(S_clu, 'clusterYPositions'), S_clu.clusterYPositions = []; end

    if isempty(S_clu.clusterXPositions) || ~isempty(S_clu.clusterYPositions)
        viClu_update = [];
    end
    if isempty(viClu_update)
        [clusterXPositions, clusterYPositions] = deal(zeros(S_clu.nClusters, 1));
        viClu1 = 1:S_clu.nClusters;
    else % selective update
        clusterXPositions = S_clu.clusterXPositions;
        clusterYPositions = S_clu.clusterYPositions;
        viClu1 = viClu_update(:)';
    end
    viSites_fet = 1:(1+P.maxSite*2-P.nSites_ref);
    for iClu = viClu1
        %     viSpk_clu1 = S_clu.spikesByCluster{iClu};
        [viSpk_clu1, viSites_clu1] = S_clu_subsample_spk_(S_clu, iClu);
        if isempty(viSpk_clu1), continue; end

        viSites_clu1 = viSites_clu1(1:end-P.nSites_ref);
        mrVp1 = squeeze_(spikeFeatures(viSites_fet,1,viSpk_clu1));
        mrSiteXY1 = single(P.mrSiteXY(viSites_clu1,:)); %electrode

        clusterXPositions(iClu) = median(centroid_mr_(mrVp1, mrSiteXY1(:,1), 2));
        clusterYPositions(iClu) = median(centroid_mr_(mrVp1, mrSiteXY1(:,2), 2));
    end
    S_clu.clusterXPositions = clusterXPositions;
    S_clu.clusterYPositions = clusterYPositions;
end %func
