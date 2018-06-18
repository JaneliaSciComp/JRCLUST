%--------------------------------------------------------------------------
function vnWav1_mean = mean_excl_(mnWav1, P)
    % calculate mean after excluding viSiteZero
    viSiteZero = get_(P, 'viSiteZero');
    if isempty(viSiteZero)
        vnWav1_mean = int16(mean(mnWav1,2));
    else
        nSites_all = size(mnWav1, 2);
        nSites_excl = numel(viSiteZero);
        nSites = nSites_all - nSites_excl;
        vnWav1_mean = int16((sum(mnWav1,2) - sum(mnWav1(:,nSites_excl),2)) / nSites);
    end
end %func
