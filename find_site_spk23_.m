%--------------------------------------------------------------------------
function [spikeSites2, spikeSites3] = find_site_spk23_(tnWav_spk, spikeSites, P)
    % find second min, excl local ref sites
    fUse_min = 0;
    imin0 = 1 - P.spkLim(1);
    viSites2 = 2:(2*P.maxSite+1-P.nSites_ref);
    miSites2 = P.miSites(viSites2,spikeSites);
    %[~,spikeSites2] = min(squeeze_(tnWav_spk(imin0,viSites2,:)));
    tnWav_spk2 = tnWav_spk(:,viSites2,:);
    if fUse_min
        mnMin_spk = squeeze_(min(tnWav_spk2));
    else
        %     mnMin_spk = -squeeze_(std(single(tnWav_spk(:,viSites2,:))));
        mnMin_spk = squeeze_(min(tnWav_spk2) - max(tnWav_spk2)); % use Vpp to determine second peak site
    end
    if nargout==1
        [~, spikeSites] = min(mnMin_spk);
        spikeSites2 = int32(mr2vr_sub2ind_(miSites2, spikeSites, []));
    else
        [~, spikePrSecSites2] = sort(mnMin_spk, 'ascend');
        spikeSites2 = int32(mr2vr_sub2ind_(miSites2, spikePrSecSites2(1,:), []));
        spikeSites3 = int32(mr2vr_sub2ind_(miSites2, spikePrSecSites2(2,:), []));
    end
end %func
