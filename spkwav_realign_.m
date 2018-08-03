%--------------------------------------------------------------------------
% 10/13/17 JJJ: Created. Realign the spikes at the min
function [spikeWaveforms1, spikeTimes1] = spkwav_realign_(spikeWaveforms1, mnWav_spk, spkLim_wav, spikeTimes1, viSite1, P)
    % subtract car and temporal shift
    % spikeWaveforms1: nSamples x nSpk x nSites_spk
    if ~strcmpi(get_set_(P, 'vcSpkRef', 'nmean'), 'nmean'), return; end

    % fprintf('\n\tRealigning spikes after LCAR (vcSpkRef=nmean)...'); t1=tic;
    dimm1 = size(spikeWaveforms1);
    trWav_spk2 = spkwav_car_(single(spikeWaveforms1), P); % apply car
    [viSpk_shift, viShift] = spkwav_shift_(trWav_spk2, 1, P);
    if isempty(viSpk_shift), return; end

    viTime_shift = spikeTimes1(viSpk_shift) - int32(viShift(:)); % spike time to shift
    spikeTimes1(viSpk_shift) = viTime_shift;
    spikeWaveforms1(:,viSpk_shift,:) = mr2tr3_(mnWav_spk, spkLim_wav, viTime_shift, viSite1);
    % fprintf('\n\t\ttook %0.1fs\n', toc(t1));
end %func
