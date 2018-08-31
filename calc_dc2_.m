%--------------------------------------------------------------------------
function dc = calc_dc2_(S0, P, vlRedo_spk)
    global spikeFeatures

    if ~all(size(spikeFeatures) ~= S0.featureDims)
        spikeFeatures = getSpikeFeatures(P);
    end

    if nargin<3, vlRedo_spk=[]; end
    fprintf('Calculating Dc\n\t'); t1=tic;
    nSites = numel(P.chanMap);
    vrDc2_site = nan(1, nSites);
    for iSite = 1:nSites
        [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(spikeFeatures, S0, P, iSite, vlRedo_spk);
        if isempty(mrFet12_), continue; end
        vrDc2_site(iSite) = compute_dc2_(mrFet12_, viiSpk12_ord_, n1_, n2_, P); % Compute DC in CPU
        fprintf('.');
    end
    dc = sqrt(abs(quantile(vrDc2_site, .5)));
    fprintf('\n\ttook %0.1fs\n', toc(t1));
end % function
