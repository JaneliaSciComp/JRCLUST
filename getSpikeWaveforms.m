%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function traces = getSpikeWaveforms(P, useRaw)
    global spikeWaveforms spikeTraces

    if nargin < 1 || isempty(P)
        P = get0_('P');
    end
    if nargin < 2
        useRaw = P.fWav_raw_show;
    end

    useRamCache = getOr(P, 'useRamCache', 1);

    if useRaw
        if ~useRamCache % clear filtered waveforms
            spikeWaveforms = [];
        end
        if isempty(spikeTraces)
            spikeTraces = load_spkraw_();
        end
        traces = spikeTraces;
    else
        if ~useRamCache % clear raw traces
            spikeTraces = [];
        end
        if isempty(spikeWaveforms)
            spikeWaveforms = load_spkwav_();
        end
        traces = spikeWaveforms;
    end
end % func
