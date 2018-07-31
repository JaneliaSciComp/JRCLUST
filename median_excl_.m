%--------------------------------------------------------------------------
function vnWav1_med = median_excl_(mnWav1, P)
    % calculate mean after excluding viSiteZero
    viSiteZero = get_(P, 'viSiteZero');
    useGPU = isGpu_(mnWav1);
    if useGPU, mnWav1 = gather_(mnWav1); end
    if isempty(viSiteZero)
        vnWav1_med = median(mnWav1,2);
    else
        viSites = setdiff(1:size(mnWav1,2), viSiteZero);
        vnWav1_med = median(mnWav1(:,viSites),2);
    end
    vnWav1_med = gpuArray_(vnWav1_med, useGPU);
end %func
