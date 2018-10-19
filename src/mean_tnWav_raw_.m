%--------------------------------------------------------------------------
function mrWav = mean_tnWav_raw_(tnWav, P)
    mrWav = meanSubt_(mean(single(tnWav),3)) * P.uV_per_bit;
end %func
