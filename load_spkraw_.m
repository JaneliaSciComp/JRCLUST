%--------------------------------------------------------------------------
% 10/10/17 JJJ: load tnWav_raw from disk
function tnWav_raw = load_spkraw_(S0)
    if nargin<1, S0 = get(0, 'UserData'); end
    tnWav_raw = load_bin_(strrep(S0.P.prmFile, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_raw);
end %func
