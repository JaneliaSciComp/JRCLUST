%--------------------------------------------------------------------------
function mrWav = mean_spikeTraces_(tnWav, P)
    mrWav = meanSubt_(mean(single(tnWav),3)) * P.uV_per_bit;
end %func
