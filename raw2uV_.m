%--------------------------------------------------------------------------
function tr = raw2uV_(spikeTraces, P)
    tr = meanSubtract(single(spikeTraces) * P.uV_per_bit);
end % function
