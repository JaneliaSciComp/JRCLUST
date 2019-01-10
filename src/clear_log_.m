%--------------------------------------------------------------------------
function S0 = clear_log_(S0)
    S0.cS_log = {};
    S0.miClu_log = [];
    set(0, 'UserData', S0);
    delete_files_(jrclust.utils.subsExt(S0.P.vcFile_prm, '_log.mat'), 0);
end %func
