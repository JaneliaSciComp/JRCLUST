%--------------------------------------------------------------------------
% 10/10/17 JJJ: load spikeTraces from disk
function spikeTraces = load_spkraw_(S0)
    if nargin<1, S0 = get(0, 'UserData'); end
    spikeTraces = load_bin_(strrep(S0.P.paramFile, '.prm', '_traces.bin'), 'int16', S0.traceDims);
end % function
