%--------------------------------------------------------------------------
function [vi1, vl1] = partition_vi_(vi, n, nBins, iBin)
    lim1 = round([iBin-1, iBin] / nBins * n) + 1;
    lim1 = min(max(lim1, 1), n+1);
    vl1 = vi >= lim1(1) & vi < lim1(2);
    vi1 = vi(vl1);
end %func
