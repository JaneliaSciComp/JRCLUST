%--------------------------------------------------------------------------
function tnWav1 = tnWav1_sites_2_(tnWav1_, viSites_spk1, viSites1, P)
    % reorder tnWav1 to viSites1

    [nT_spk, nSites_spk, nSpk1] = size(tnWav1_);
    % nSites1 = numel(viSites1);
    nSites = numel(P.viSite2Chan);
    tnWav1 = zeros([nT_spk, nSites, nSpk1], 'like', tnWav1_); %full
    % assert_(nSites_spk==nSites1, 'tnWav1_sites_: nSites must agree');
    % tnWav1 = zeros([nT_spk, nSpk1, nSites1], 'like', tnWav1_); %subset of spk, complete
    for iSite1 = 1:numel(viSites1)
        iSite11 = viSites1(iSite1);
        viSpk1 = find(viSites_spk1 == iSite11);
        if isempty(viSpk1), return; end
        viSites11 = P.miSites(:, iSite11);
        tnWav1(:,viSites11,viSpk1) = tnWav1_(:,:,viSpk1);
    end
    % tnWav1 = permute(tnWav1, [1,3,2]);
    tnWav1 = tnWav1(:, viSites1, :); %reduced subset
end %func
