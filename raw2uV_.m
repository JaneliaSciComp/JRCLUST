%--------------------------------------------------------------------------
function tr = raw2uV_(spikeTraces, P)
    tr = meanSubt_(single(spikeTraces) * P.uV_per_bit);
end %func
