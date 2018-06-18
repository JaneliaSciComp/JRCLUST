%--------------------------------------------------------------------------
function mr_ref = tr_sort_ref_(trWav_spk1, viSites_ref)
    %  trWav_spk1: nT x nSpk x nSites, single
    fUseMin = 0;
    dimm1 = size(trWav_spk1);
    tr0 = gather_(trWav_spk1);
    %tr = permute(gather_(trWav_spk1), [3,1,2]);
    mr_ref = zeros(dimm1([1,2]), 'like', tr0);
    if fUseMin
        P = get0_('P');
        iT0 = 1 - P.spkLim(1);
        [~, miSites_ref] = sort(permute(tr0(iT0,:,:), [3,2,1]), 'ascend');
    else % use std
        [~, miSites_ref] = sort(permute(var(tr0), [3,2,1]), 'descend'); % use lest activities for ref
    end
    miSites_ref = miSites_ref(viSites_ref,:);
    for iSpk1=1:dimm1(2)
        %mr_ref(:,iSpk1) = mean(tr(miSites_ref(:,iSpk1),:,iSpk1));
        mr_ref(:,iSpk1) = mean(tr0(:,iSpk1,miSites_ref(:,iSpk1)), 3);
    end
end %func
