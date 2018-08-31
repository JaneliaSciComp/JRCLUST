%--------------------------------------------------------------------------
% 10/10/17 JJJ: load spikeWaveforms from disk
function spikeWaveforms = load_spkwav_(S0)
    if nargin<1, S0 = get(0, 'UserData'); end
    spikeWaveforms = load_bin_(strrep(S0.P.paramFile, '.prm', '_waveforms.bin'), 'int16', S0.waveformDims);
end % function
