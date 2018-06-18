%--------------------------------------------------------------------------
function vr = quantile_(mr, p)
    n = size(mr,1);
    idx = max(min(round(n * p), n), 1);
    mr = sort(mr, 'ascend');
    vr = mr(idx,:);
end %func
