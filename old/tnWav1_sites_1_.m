%--------------------------------------------------------------------------
function tnWav1 = tnWav1_sites_1_(tnWav1_, miSites1, viSites1)
    [nT_spk, nSites_spk, nSpk1] = size(tnWav1_);
    nSites1 = numel(viSites1);
    % assert_(nSites_spk==nSites1, 'tnWav1_sites_: nSites must agree');
    tnWav1 = zeros([nT_spk, nSpk1, nSites1], 'like', tnWav1_); %subset of spk, complete
    for iSite1 = 1:nSites1
        [viSite11, viiSpk11] = find(miSites1 == viSites1(iSite1));
        nSpk11 = numel(viiSpk11);
        mnWav_spk11 = reshape(tnWav1_(:, :, viiSpk11), nT_spk, []);
        mnWav_spk11 = mnWav_spk11(:, sub2ind([nSites_spk, nSpk11], viSite11', 1:nSpk11));
        tnWav1(:, viiSpk11, iSite1) = mnWav_spk11;
    end
    tnWav1 = permute(tnWav1, [1,3,2]);
end %func
