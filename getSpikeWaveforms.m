%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function traces = getSpikeWaveforms(P, useRaw)
    % get (raw or filtered) spike waveforms
    % return value: nSamples x nSites x nSpikes tensor

    global spikeWaveforms spikeTraces

    S0 = get(0, 'UserData');

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
        if ~allDimsEqual(spikeTraces, S0.traceDims)
            spikeTraces = load_spkraw_();
        end
        traces = spikeTraces;
    else
        if ~useRamCache % clear raw traces
            spikeTraces = [];
        end
        if ~allDimsEqual(spikeWaveforms, S0.waveformDims)
            spikeWaveforms = load_spkwav_();
        end
        traces = spikeWaveforms;
    end
end % function
