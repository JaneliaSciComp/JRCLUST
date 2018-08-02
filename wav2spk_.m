%--------------------------------------------------------------------------
function [tnWav_spk_raw, tnWav_spk, trFet_spk, miSite_spk, spikeTimes, vnAmp_spk, siteThresholds, useGPU] = ...
    wav2spk_(mnWav1, vrWav_mean1, P, spikeTimes, spikeSites, mnWav1_pre, mnWav1_post)
    % tnWav_spk: spike waveform. nSamples x nSites x nSpikes
    % trFet_spk: nSites x nSpk x nFet
    % miSite_spk: nSpk x nFet
    % spikes are ordered in time
    % spikeSites and spikeTimes is uint32 format, and tnWav_spk: single format
    % mnWav1: raw waveform (unfiltered)
    % wav2spk_(mnWav1, vrWav_mean1, P)
    % wav2spk_(mnWav1, vrWav_mean1, P, spikeTimes, spikeSites)
    % 6/27/17 JJJ: accurate spike detection at the overlap region
    % 6/29/17 JJJ: matched filter supported

    if nargin<4, spikeTimes = []; end
    if nargin<5, spikeSites = []; end
    if nargin<6, mnWav1_pre = []; end
    if nargin<7, mnWav1_post = []; end
    [tnWav_spk_raw, tnWav_spk, trFet_spk, miSite_spk] = deal([]);
    nFet_use = get_set_(P, 'nFet_use', 2);
    fMerge_spk = 1; %debug purpose
    fShift_pos = 0; % shift center position based on center of mass
    % fRecenter_spk = 0;
    nSite_use = P.maxSite*2+1 - P.nSites_ref;
    if nSite_use==1, nFet_use=1; end
    siteThresholds = get0_('siteThresholds');
    nPad_pre = size(mnWav1_pre,1);

    %-----
    % Filter
    fprintf('\tFiltering spikes...'); t_filter = tic;
    if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)
        mnWav1 = [mnWav1_pre; mnWav1; mnWav1_post];
    end
    % [mnWav2, vnWav11, mnWav1, P.useGPU] = wav_preproces_(mnWav1, P);
    mnWav1_ = mnWav1; % keep a copy in CPU
    try
        [mnWav1, P.useGPU] = gpuArray_(mnWav1, P.useGPU);
        if P.fft_thresh>0, mnWav1 = fft_clean_(mnWav1, P); end
        [mnWav2, vnWav11] = filt_car_(mnWav1, P);
    catch % GPu failure
        P.useGPU = 0;
        mnWav1 = mnWav1_;
        if P.fft_thresh>0, mnWav1 = fft_clean_(mnWav1, P); end
        [mnWav2, vnWav11] = filt_car_(mnWav1, P);
    end
    mnWav1_ = []; %remove from memory


    %-----
    % common mode rejection
    if P.blank_thresh > 0
        if isempty(vnWav11)
            vnWav11 = mr2ref_(mnWav2, P.vcCommonRef, P.viSiteZero); %vrWav_mean1(:);
        end
        vlKeep_ref = car_reject_(vnWav11(:), P);
        fprintf('Rejecting %0.3f %% of time due to motion\n', (1-mean(vlKeep_ref))*100 );
    else
        vlKeep_ref = [];
    end
    % setUserData(vlKeep_ref);
    fprintf('\ttook %0.1fs\n', toc(t_filter));

    switch get_set_(P, 'vcFilter_detect', '')
        case {'', 'none'}, mnWav3 = mnWav2;
        case 'ndist'
        [mnWav3, nShift_post] = filter_detect_(mnWav1, P); % pass raw trace
        otherwise
        [mnWav3, nShift_post] = filter_detect_(mnWav2, P); % pass filtered trace
    end

    %-----
    % detect spikes or use the one passed from the input (importing)
    if isempty(siteThresholds)
        try
            siteThresholds = gather_(int16(mr2rms_(mnWav3, 1e5) * P.qqFactor));
        catch
            siteThresholds = int16(mr2rms_(gather_(mnWav3), 1e5) * P.qqFactor);
            P.useGPU = 0;
        end
    end
    if isempty(spikeTimes) || isempty(spikeSites)
        P_ = setfield(P, 'nPad_pre', nPad_pre);
        [spikeTimes, vnAmp_spk, spikeSites] = detect_spikes_(mnWav3, siteThresholds, vlKeep_ref, P_);
    else
        spikeTimes = spikeTimes + nPad_pre;
        vnAmp_spk = mnWav3(sub2ind(size(mnWav3), spikeTimes, spikeSites)); % @TODO read spikes at the site and time
    end
    vnAmp_spk = gather_(vnAmp_spk);
    % if nShift_post~=0, spikeTimes = spikeTimes + nShift_post; end % apply possible shift due to filtering

    % reject spikes within the overlap region
    if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)
        ilim_spk = [nPad_pre+1, size(mnWav3,1) - size(mnWav1_post,1)]; %inclusive
        viKeep_spk = find(spikeTimes >= ilim_spk(1) & spikeTimes <= ilim_spk(2));
        [spikeTimes, vnAmp_spk, spikeSites] = multifun_(@(x)x(viKeep_spk), spikeTimes, vnAmp_spk, spikeSites);
    end %if
    if isempty(spikeTimes), return; end


    %-----
    % Extract spike waveforms and build a spike table
    fprintf('\tExtracting features'); t_fet = tic;
    % mnWav2 = gather_(mnWav2); %do in CPU. 10.2s in GPU, 10.4s in CPU
    % if fRecenter_spk % center site is where the energy is the highest, if disabled min is chosen
    %     tnWav_spk = mn2tn_wav_spk2_(mnWav2, spikeSites, spikeTimes, P);
    %     %[~, viMaxSite_spk] = max(squeeze_(std(single(tnWav_spk))));
    %     [~, viMaxSite_spk] = max(squeeze_(max(tnWav_spk) - min(tnWav_spk)));
    %     spikeSites = P.miSites(sub2ind(size(P.miSites), viMaxSite_spk(:), spikeSites));
    % end
    spikeSites_ = gpuArray_(spikeSites);
    [tnWav_spk_raw, tnWav_spk, spikeTimes] = mn2tn_wav_(mnWav1, mnWav2, spikeSites_, spikeTimes, P); fprintf('.');
    if nFet_use >= 2
        viSite2_spk = find_site_spk23_(tnWav_spk, spikeSites_, P);
        tnWav_spk2 = mn2tn_wav_spk2_(mnWav2, viSite2_spk, spikeTimes, P);
    else
        [viSite2_spk, tnWav_spk2] = deal([]);
    end

    %-----
    % Cancel overlap
    if get_set_(P, 'fCancel_overlap', 0)
        try
            [tnWav_spk, tnWav_spk2] = cancel_overlap_spk_(tnWav_spk, tnWav_spk2, spikeTimes, spikeSites, viSite2_spk, siteThresholds, P);
        catch
            fprintf(2, 'fCancel_overlap failed\n');
        end
    end

    tnWav_spk_raw = gather_(tnWav_spk_raw);
    dialogAssert(nSite_use >0, 'nSites_use = maxSite*2+1 - nSites_ref must be greater than 0');
    switch nFet_use
        case 3
        [viSite2_spk, viSite3_spk] = find_site_spk23_(tnWav_spk, spikeSites_, P); fprintf('.');
        mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
        mrFet2 = trWav2fet_(tnWav_spk2, P); fprintf('.');
        mrFet3 = trWav2fet_(mn2tn_wav_spk2_(mnWav2, viSite3_spk, spikeTimes, P), P); fprintf('.');
        trFet_spk = permute(cat(3, mrFet1, mrFet2, mrFet3), [1,3,2]); %nSite x nFet x nSpk
        miSite_spk = [spikeSites_(:), viSite2_spk(:), viSite3_spk(:)]; %nSpk x nFet
        case 2
        mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
        mrFet2 = trWav2fet_(tnWav_spk2, P); fprintf('.');
        trFet_spk = permute(cat(3, mrFet1, mrFet2), [1,3,2]); %nSite x nFet x nSpk
        miSite_spk = [spikeSites_(:), viSite2_spk(:)]; %nSpk x nFet
        case 1
        mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
        trFet_spk = permute(mrFet1, [1,3,2]); %nSite x nFet x nSpk
        miSite_spk = [spikeSites_(:)];
        otherwise
        error('wav2spk_: nFet_use must be 1, 2 or 3');
    end

    if nPad_pre > 0, spikeTimes = spikeTimes - nPad_pre; end
    [spikeTimes, trFet_spk, miSite_spk, tnWav_spk] = ...
    gather_(spikeTimes, trFet_spk, miSite_spk, tnWav_spk);
    useGPU = P.useGPU;
    fprintf('\ttook %0.1fs\n', toc(t_fet));
end %func
