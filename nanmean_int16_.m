%--------------------------------------------------------------------------
function mrWav_clu1 = nanmean_int16_(tnWav0, dimm_mean, fUseCenterSpk, iSite1, viSite0, P); % * S0.P.uV_per_bit;
    if fUseCenterSpk
        mrWav_clu1 = mean(single(tnWav0), dimm_mean);
    else
        %     nSites = numel(P.viSite2Chan);
        viSite1 = P.miSites(:, iSite1);
        trWav = nan([size(tnWav0,1), numel(viSite1), numel(viSite0)], 'single');
        viSites_uniq = unique(viSite0);
        nSites_uniq = numel(viSites_uniq);
        miSites_uniq = P.miSites(:, viSites_uniq);
        for iSite_uniq1 = 1:nSites_uniq
            iSite_uniq = viSites_uniq(iSite_uniq1);
            viSpk_ = find(viSite0 == iSite_uniq);
            [~, viSite1a_, viSite1b_] = intersect(viSite1, miSites_uniq(:,iSite_uniq1));
            if isempty(viSite1a_), continue; end
            trWav(:, viSite1a_, viSpk_) = tnWav0(:,viSite1b_,viSpk_); % P.miSites(:,iSite_unique)
        end
        mrWav_clu1 = nanmean(trWav, dimm_mean);
        %     mrWav_clu1 = nanmean(trWav(:,P.miSites(:, iSite1),:), dimm_mean);
    end
    mrWav_clu1 = meanSubt_(mrWav_clu1); %122717 JJJ
end %func
