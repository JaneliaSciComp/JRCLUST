%--------------------------------------------------------------------------
function min_ = min_step_(vr)
    min_ = diff(sort(vr));
    min_ = min(min_(min_>0));
end %func
