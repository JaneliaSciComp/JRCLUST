%--------------------------------------------------------------------------
% 10/10/17 JJJ: load tnWav_spk from disk
function tnWav_spk = load_spkwav_(S0)
    if nargin<1, S0 = get(0, 'UserData'); end
    tnWav_spk = load_bin_(strrep(S0.P.vcFile_prm, '.prm', '_spkwav.jrc'), 'int16', S0.dimm_spk);
end %func
