%--------------------------------------------------------------------------
function [tnWav_spk, vnThresh_site] = wav2spk_gt_(mnWav1, P, viTime_spk, mnWav1_pre, mnWav1_post)
    % tnWav_spk: spike waveform. nSamples x nSites x nSpikes
    % trFet_spk: nSites x nSpk x nFet
    % miSite_spk: nSpk x nFet
    % spikes are ordered in time
    % viSite_spk and viTime_spk is uint32 format, and tnWav_spk: single format
    % mnWav1: raw waveform (unfiltered)
    % wav2spk_(mnWav1, vrWav_mean1, P)
    % wav2spk_(mnWav1, vrWav_mean1, P, viTime_spk, viSite_spk)
    % 7/5/17 JJJ: accurate spike detection at the overlap region

    if nargin<5, mnWav1_pre = []; end
    if nargin<6, mnWav1_post = []; end

    % Filter
    if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)
        mnWav1 = [mnWav1_pre; mnWav1; mnWav1_post];
    end
    if P.fft_thresh>0, mnWav1 = jrclust.filters.fftClean(mnWav1, P); end
    [mnWav2, vnWav11] = jrclust.filters.filtCAR(mnWav1, P); % filter and car

    % detect spikes or use the one passed from the input (importing)
    vnThresh_site = jrclust.utils.tryGather(int16(jrclust.utils.estimateRMS(mnWav2, 1e5) * P.qqFactor));
    nPad_pre = size(mnWav1_pre,1);
    viTime_spk = viTime_spk + nPad_pre;

    % reject spikes within the overlap region: problem if viTime prespecified.
    if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)
        ilim_spk = [nPad_pre+1, size(mnWav2,1) - size(mnWav1_post,1)]; %inclusive
        viKeep_spk = find(viTime_spk >= ilim_spk(1) & ilim_spk <= ilim_spk(2));
        viTime_spk = multifun_(@(x)x(viKeep_spk), viTime_spk);
    end %if
    % [tnWav_spk_raw, tnWav_spk] = mn2tn_wav_(mnWav1, mnWav2, [], viTime_spk, P);
    tnWav_spk = permute(jrclust.utils.tryGather(jrclust.utils.extractWindows(mnWav2, P.spkLim, viTime_spk)), [1,3,2]);

    % if nPad_pre > 0, viTime_spk = viTime_spk - nPad_pre; end
end %func
