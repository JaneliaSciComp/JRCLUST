%--------------------------------------------------------------------------
function [spikeSites2, spikeSites3] = find_site_spk23_(spikeWaveforms, spikeSites, P)
    % find second min, excl local ref sites
    fUse_min = 0;
    imin0 = 1 - P.spkLim(1);
    viSites2 = 2:(2*P.maxSite+1-P.nSites_ref);
    miSites2 = P.miSites(viSites2,spikeSites);
    %[~,spikeSites2] = min(squeeze_(spikeWaveforms(imin0,viSites2,:)));
    spikeWaveforms2 = spikeWaveforms(:,viSites2,:);
    if fUse_min
        mnMin_spk = squeeze_(min(spikeWaveforms2));
    else
        %     mnMin_spk = -squeeze_(std(single(spikeWaveforms(:,viSites2,:))));
        mnMin_spk = squeeze_(min(spikeWaveforms2) - max(spikeWaveforms2)); % use Vpp to determine second peak site
    end
    if nargout==1
        [~, spikeSites] = min(mnMin_spk);
        spikeSites2 = int32(mr2vr_sub2ind_(miSites2, spikeSites, []));
    else
        [~, spikePrimarySecondarySites2] = sort(mnMin_spk, 'ascend');
        spikeSites2 = int32(mr2vr_sub2ind_(miSites2, spikePrimarySecondarySites2(1,:), []));
        spikeSites3 = int32(mr2vr_sub2ind_(miSites2, spikePrimarySecondarySites2(2,:), []));
    end
end % function
