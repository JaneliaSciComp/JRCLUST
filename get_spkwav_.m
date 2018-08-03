%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function tnWav_ = get_spkwav_(P, fRaw)
    % if ~fRamCache, only keep one of the spikeTraces or spikeWaveforms in memory
    global spikeWaveforms spikeTraces
    if nargin<1, P = []; end
    if isempty(P), P = get0_('P'); end
    if nargin<2, fRaw = P.fWav_raw_show; end

    fRamCache = get_set_(P, 'fRamCache', 1);
    if fRaw
        if ~fRamCache, spikeWaveforms = []; end % clear spk
        if isempty(spikeTraces), spikeTraces = load_spkraw_(); end
        tnWav_ = spikeTraces;
    else
        if ~fRamCache, spikeTraces = []; end % clear raw
        if isempty(spikeWaveforms), spikeWaveforms = load_spkwav_(); end
        tnWav_ = spikeWaveforms;
    end
end %func
