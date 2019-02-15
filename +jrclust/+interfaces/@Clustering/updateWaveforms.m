function updateWaveforms(obj, updateMe)
    %UPDATEWAVEFORMS Update mean waveforms and sim scores
    if nargin < 2
        updateMe = [];
    end

    obj.computeMeanWaveforms(updateMe);
    obj.computeWaveformSim(updateMe);
end