%--------------------------------------------------------------------------
function vi = subsample_vr_(vi, nMax)
    if numel(vi)>nMax
        nSkip = floor(numel(vi)/nMax);
        if nSkip>1, vi = vi(1:nSkip:end); end
        if numel(vi)>nMax, vi = vi(1:nMax); end
    end
end %func
