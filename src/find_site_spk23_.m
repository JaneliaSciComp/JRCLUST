%--------------------------------------------------------------------------
function [viSite_spk2, viSite_spk3] = find_site_spk23_(tnWav_spk, viSite_spk, P)
    % find second min, excl local ref sites
    fUse_min = 0;
    imin0 = 1 - P.spkLim(1);
    viSites2 = 2:(2*P.maxSite+1-P.nSites_ref);
    miSites2 = P.miSites(viSites2,viSite_spk);
    %[~,viSite_spk2] = min(squeeze_(tnWav_spk(imin0,viSites2,:)));
    tnWav_spk2 = tnWav_spk(:,viSites2,:);
    if fUse_min
        mnMin_spk = squeeze_(min(tnWav_spk2));
    else
        %     mnMin_spk = -squeeze_(std(single(tnWav_spk(:,viSites2,:))));
        mnMin_spk = squeeze_(min(tnWav_spk2) - max(tnWav_spk2)); % use Vpp to determine second peak site
    end
    if nargout==1
        if numel(viSites2) == 1
            viSite_spk2 = int32(miSites2);
        else
            [~, viSite_spk] = min(mnMin_spk, [], 1);
            viSite_spk2 = int32(jrclust.utils.rowColSelect(miSites2, viSite_spk, []));
        end
    else
        [~, miSite_spk2] = sort(mnMin_spk, 'ascend');
        viSite_spk2 = int32(jrclust.utils.rowColSelect(miSites2, miSite_spk2(1,:), []));
        viSite_spk3 = int32(jrclust.utils.rowColSelect(miSites2, miSite_spk2(2,:), []));
    end
end %func
