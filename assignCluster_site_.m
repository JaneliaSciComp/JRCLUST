%--------------------------------------------------------------------------
function viClu = assignCluster_site_(S_clu, S0)
    nRepeat = 100;

    nSites = numel(S0.cviSpk_site);
    viClu = zeros(numel(S_clu.rho), 1, 'int32');
    viClu(S_clu.clusterCenters) = 1:numel(S_clu.clusterCenters);
    % for iRepeat = 1:nRepeat
    for iSite = 1:nSites
        viSpk_ = S0.cviSpk_site{iSite}; %find(S0.spikeSites == iSite);
        if isempty(viSpk_), continue; end
        %     viSpk_ = viSpk_(S0.spikeSites(S_clu.nneigh(viSpk_)) == iSite); % in group spikes only
        %     viSpk_ = viSpk_(ismember(S_clu.nneigh(viSpk_), viSpk_));
        cl_ = viClu(viSpk_);
        ordrho_ = rankorder_(S_clu.rho(viSpk_), 'descend');
        [vl_, nneigh_] = ismember(S_clu.nneigh(viSpk_), viSpk_);
        vi_ = find(vl_);
        vi_ = vi_(:)';
        for iRepeat = 1:nRepeat
            vi_(cl_(ordrho_(vi_))~=0) = [];
            if isempty(vi_), break; end
            for i_ = vi_
                cl_(ordrho_(i_)) = cl_(nneigh_(i_));
            end
        end
        viClu(viSpk_) = cl_;
    end %for
    %     fprintf('%d: %f\n', iRepeat, mean(viClu>0));
end % function
