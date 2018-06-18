%--------------------------------------------------------------------------
function trFet_spk_ = nneigh_ave_(S_clu, P, trFet_spk)
    % cluster based averaging fet cleanup
    error('not done');
    nClu = S_clu.nClu;
    S0 = get(0, 'UserData');
    trFet_spk = gpuArray(trFet_spk);
    trFet_spk_ = zeros(size(trFet_spk), 'like', trFet_spk);
    for iClu = 1:nClu
        viSpk1 = S_clu.cviSpk_clu{iClu};
        trFet1 = permute(trFet_spk(:,:,viSpk1), [1,3,2]);
        viSite1 = S0.viSite_spk(viSpk1);
        viSite2 = S0.viSite2_spk(viSpk1);
        viSite_unique1 = unique(viSite1);
        for iSite1 = 1:numel(viSite_unique1)
            iSite11 = viSite_unique1(iSite1);
            viSpk11 = find(viSite1 == iSite11);
            viSpk12 = find(viSite2 == iSite11);
            mrFet11 = trFet1(:,viSpk11,1);
            mrFet12 = [mrFet11, trFet1(:,viSpk12,2)];
            [vr11, vi_nneigh12] = min(set_diag_(eucl2_dist_(mrFet12, mrFet11), nan));
            %         trFet_spk_(:,:,
        end
    end
end %func
