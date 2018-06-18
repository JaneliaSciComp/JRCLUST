%--------------------------------------------------------------------------
function vrSnr_clu = S_clu_snr_(S_clu)
    S0 = get(0, 'UserData');
    S_clu = S_clu_wav_(S_clu, [], 1);
    mrVmin_clu = shiftdim(min(S_clu.tmrWav_clu,[],1));
    [vrVmin_clu, viSite_clu] = min(mrVmin_clu,[],1);
    vrVrms_site = single(S0.vrThresh_site(:)) / S0.P.qqFactor;
    vrSnr_clu = abs(vrVmin_clu(:)) ./ bit2uV_(vrVrms_site(viSite_clu), S0.P);
end %func
