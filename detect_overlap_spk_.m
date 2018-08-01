%--------------------------------------------------------------------------
% 12/16/17 JJJ: find overlapping spikes. only return spikes more negative than others
function [viSpk_ol_spk, vnDelay_ol_spk, vnCount_ol_spk] = detect_overlap_spk_(spikeTimes, viSite_spk, P);

    mrDist_site = squareform(pdist(P.mrSiteXY));
    nlimit = int32(diff(P.spkLim));
    maxDist_site_um = P.maxDist_site_um;
    nSpk = numel(spikeTimes);
    nSites = max(viSite_spk);
    cviSpk_site = arrayfun(@(iSite)int32(find(viSite_spk==iSite)), (1:nSites)', 'UniformOutput', 0);
    spikeTimes = gather_(spikeTimes);
    [viSpk_ol_spk, vnDelay_ol_spk, vnCount_ol_spk] = deal(zeros(size(viSite_spk), 'int32'));
    for iSite = 1:nSites
        viSpk1 = cviSpk_site{iSite};
        if isempty(viSpk1), continue; end
        viSite2 = find(mrDist_site(:,iSite) <= maxDist_site_um & mrDist_site(:,iSite) > 0);
        viSpk2 = cell2mat_(cviSpk_site(viSite2));
        [n1, n2] = deal(numel(viSpk1), numel(viSpk2));
        viSpk12 = [viSpk1(:); viSpk2(:)];
        [viTime1, viTime12] = deal(spikeTimes(viSpk1), spikeTimes(viSpk12));

        % find overlapping spikes that has smaller amplitudes and within site limit
        %     [viOverlap1, viDelay1] = deal(zeros(size(viTime12), 'like', viTime12));
        %     vlOverlap1 = false(size(viTime1));
        for iDelay = 0:nlimit
            [vl12_, vi1_] = ismember(viTime12, viTime1 + iDelay);
            if iDelay == 0 , vl12_(1:n1) = 0; end % exclude same site comparison
            vi12_ = find(vl12_);
            if isempty(vi12_), continue; end
            [viSpk1_, viSpk12_] = deal(viSpk1(vi1_(vi12_)), viSpk12(vi12_));
            %         viiSpk_ = find(spikeTimes(viSpk12_) < spikeTimes(viSpk1_)); % pick earlier time only
            %         if isempty(viiSpk_), continue; end;
            %         [viSpk1_, viSpk12_] = deal(viSpk1_(viiSpk_), viSpk12_(viiSpk_));
            viSpk_ol_spk(viSpk12_) = viSpk1_;
            vnDelay_ol_spk(viSpk12_) = iDelay;
            vnCount_ol_spk(viSpk12_) = vnCount_ol_spk(viSpk12_) + 1; % 13% spikes collide
        end
    end
end %func
