%--------------------------------------------------------------------------
function [vr1, vi] = subsample_(vr, nSubsample)
    vi = subsample_vr_(1:numel(vr), nSubsample);
    vr1 = vr(vi);
end % function
