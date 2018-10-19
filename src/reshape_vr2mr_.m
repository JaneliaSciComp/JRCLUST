%--------------------------------------------------------------------------
function mr = reshape_vr2mr_(vr, nwin)
    nbins = ceil(numel(vr)/nwin);
    vr(nbins*nwin) = 0; %expand size
    mr = reshape(vr(1:nbins*nwin), nwin, nbins);
end %func
