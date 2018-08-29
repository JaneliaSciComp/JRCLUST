%--------------------------------------------------------------------------
function [spikeWaveforms, siteThresholds] = wav2spk_gt_(mnWav1, P, spikeTimes, mnWav1_pre, mnWav1_post)
    % spikeWaveforms: spike waveform. nSamples x nSites x nSpikes
    % spikeFeatures: nSites x nSpk x nFet
    % spikePrimarySecondarySites: nSpk x nFet
    % spikes are ordered in time
    % spikeSites and spikeTimes is uint32 format, and spikeWaveforms: single format
    % mnWav1: raw waveform (unfiltered)
    % wav2spk_(mnWav1, vrWav_mean1, P)
    % wav2spk_(mnWav1, vrWav_mean1, P, spikeTimes, spikeSites)
    % 7/5/17 JJJ: accurate spike detection at the overlap region

    if nargin<5, mnWav1_pre = []; end
    if nargin<6, mnWav1_post = []; end

    % Filter
    if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)
        mnWav1 = [mnWav1_pre; mnWav1; mnWav1_post];
    end
    if P.fft_thresh>0, mnWav1 = fft_clean_(mnWav1, P); end
    [mnWav2, vnWav11] = filt_car_(mnWav1, P); % filter and car

    % detect spikes or use the one passed from the input (importing)
    siteThresholds = gather_(int16(mr2rms_(mnWav2, 1e5) * P.qqFactor));
    nPad_pre = size(mnWav1_pre,1);
    spikeTimes = spikeTimes + nPad_pre;

    % reject spikes within the overlap region: problem if viTime prespecified.
    if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)
        ilim_spk = [nPad_pre+1, size(mnWav2,1) - size(mnWav1_post,1)]; %inclusive
        viKeep_spk = find(spikeTimes >= ilim_spk(1) & ilim_spk <= ilim_spk(2));
        spikeTimes = multifun_(@(x)x(viKeep_spk), spikeTimes);
    end %if
    % [spikeWaveforms_raw, spikeWaveforms] = mn2tn_wav_(mnWav1, mnWav2, [], spikeTimes, P);
    spikeWaveforms = permute(gather_(mr2tr3_(mnWav2, P.spkLim, spikeTimes)), [1,3,2]);

    % if nPad_pre > 0, spikeTimes = spikeTimes - nPad_pre; end
end % function
