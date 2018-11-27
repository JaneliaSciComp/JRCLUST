%--------------------------------------------------------------------------
function tr = raw2uV_(tnWav_raw, P)
    tr = jrclust.utils.meanSubtract(single(tnWav_raw) * P.uV_per_bit);
end %func
