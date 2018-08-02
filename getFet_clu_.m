%--------------------------------------------------------------------------
function [mrFet1, viSpk1] = getFet_clu_(iClu1, iSite, S0)
    % get features on-fly
    MAX_SAMPLE = 10000; % max points ti display
    if nargin<3
        [S_clu, P, viSite_spk] = get0_('S_clu', 'P', 'viSite_spk');
    else
        [S_clu, P, viSite_spk] = deal(S0.S_clu, S0.P, S0.viSite_spk);
    end
    % if nargin<2, viSite = P.miSites(:, S0.S_clu.clusterSites(iClu1)); end
    if isempty(iClu1) % select spikes based on sites
        n_use = 1 + round(P.maxSite);
        viSite_ = P.miSites(1:n_use, iSite);
        try
            viSpk1 = cell2mat_([S0.cviSpk_site(viSite_)]');
        catch
            viSpk1 = find(ismember(viSite_spk, viSite_));
        end
        viSpk1 = randomSelect_(viSpk1, MAX_SAMPLE);
    else
        viSpk1 = S_clu.cviSpk_clu{iClu1};
    end

    switch lower(P.vcFet_show)
        case {'vmin', 'vpp'}
        mrWav_spk1 = squeeze_(tnWav2uV_(tnWav_spk_sites_(viSpk1, iSite, S0), P));
        mrFet1 = max(mrWav_spk1)-min(mrWav_spk1);
        case 'cov'
        mrFet1 = calc_cov_spk_(viSpk1, iSite);
        case {'pca', 'gpca'}
        mrFet1 = pca_pc_spk_(viSpk1, iSite);
        case {'ppca', 'private pca'}
        [mrPv1, mrPv2] = pca_pv_clu_(iSite, S0.iCluCopy);
        mrFet1 = pca_pc_spk_(viSpk1, iSite, mrPv1, mrPv2);
        otherwise
        error('not implemented yet');
    end
    mrFet1 = squeeze_(abs(mrFet1));
end %func
