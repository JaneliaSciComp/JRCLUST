%--------------------------------------------------------------------------
function [spikeTraces, spikeWaveforms, spikeFeatures, spikePrimarySecondarySites, spikeTimes, vnAmp_spk, siteThresholds, useGPU] = ...
    wav2spk_(rawTraces, channelMeans, P, spikeTimes, spikeSites, prePadding, postPadding)
    % spikeWaveforms: spike waveform. nSamples x nSites x nSpikes
    % spikeFeatures: nSites x nSpk x nFet
    % spikePrimarySecondarySites: nSpk x nFet
    % spikes are ordered in time
    % spikeSites and spikeTimes is uint32 format, and spikeWaveforms: single format
    % rawTraces: raw waveform (unfiltered)
    % wav2spk_(rawTraces, channelMeans, P)
    % wav2spk_(rawTraces, channelMeans, P, spikeTimes, spikeSites)
    % 6/27/17 JJJ: accurate spike detection at the overlap region
    % 6/29/17 JJJ: matched filter supported

    if nargin < 4
        spikeTimes = [];
    end
    if nargin < 5
        spikeSites = [];
    end
    if nargin < 6
        prePadding = [];
    end
    if nargin < 7
        postPadding = [];
    end

    [spikeTraces, spikeWaveforms, spikeFeatures, spikePrimarySecondarySites] = deal([]);
    nFet_use = getOr(P, 'nFet_use', 2);

    nSite_use = P.maxSite*2+1 - P.nSites_ref;
    if nSite_use == 1
        nFet_use = 1;
    end

    siteThresholds = get0_('siteThresholds');
    nPad_pre = size(prePadding, 1);

    %-----
    % Filter
    fprintf('\tFiltering spikes...');
    t_filter = tic;

    rawTraces = [prePadding; rawTraces; postPadding];
    rawTracesCPU = rawTraces; % keep a copy in CPU
    try
        [rawTraces, P.useGPU] = gpuArray_(rawTraces, P.useGPU);
    catch % GPu failure
        P.useGPU = 0;
        rawTraces = rawTracesCPU;
    end

    if isfield(P, 'fft_thresh') && P.fft_thresh > 0
        rawTraces = fft_clean_(rawTraces, P);
    end

    [filteredTraces, vnWav11] = filt_car_(rawTraces, P);
    rawTracesCPU = []; % clear from memory

    %-----
    % common mode rejection
    if P.blank_thresh > 0
        if isempty(vnWav11)
            vnWav11 = mr2ref_(filteredTraces, P.vcCommonRef, P.viSiteZero); %channelMeans(:);
        end

        vlKeep_ref = car_reject_(vnWav11(:), P);
        fprintf('Rejecting %0.3f %% of time due to motion\n', (1-mean(vlKeep_ref))*100 );
    else
        vlKeep_ref = [];
    end
    % setUserData(vlKeep_ref);
    fprintf('\ttook %0.1fs\n', toc(t_filter));

    switch getOr(P, 'vcFilter_detect', '')
        case {'', 'none'}
            mnWav3 = filteredTraces;
        case 'ndist'
            [mnWav3, nShift_post] = filter_detect_(rawTraces, P); % pass raw trace
        otherwise
            [mnWav3, nShift_post] = filter_detect_(filteredTraces, P); % pass filtered trace
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
    if ~isempty(prePadding) || ~isempty(postPadding)
        ilim_spk = [nPad_pre+1, size(mnWav3,1) - size(postPadding,1)]; %inclusive
        viKeep_spk = find(spikeTimes >= ilim_spk(1) & spikeTimes <= ilim_spk(2));
        [spikeTimes, vnAmp_spk, spikeSites] = multifun_(@(x) x (viKeep_spk), spikeTimes, vnAmp_spk, spikeSites);
    end % if

    if isempty(spikeTimes)
        return;
    end


    %-----
    % Extract spike waveforms and build a spike table
    fprintf('\tExtracting features'); t_fet = tic;
    % filteredTraces = gather_(filteredTraces); %do in CPU. 10.2s in GPU, 10.4s in CPU
    % if fRecenter_spk % center site is where the energy is the highest, if disabled min is chosen
    %     spikeWaveforms = mn2tn_wav_spk2_(filteredTraces, spikeSites, spikeTimes, P);
    %     %[~, viMaxSite_spk] = max(squeeze_(std(single(spikeWaveforms))));
    %     [~, viMaxSite_spk] = max(squeeze_(max(spikeWaveforms) - min(spikeWaveforms)));
    %     spikeSites = P.miSites(sub2ind(size(P.miSites), viMaxSite_spk(:), spikeSites));
    % end
    spikeSites_ = gpuArray_(spikeSites);
    [spikeTraces, spikeWaveforms, spikeTimes] = mn2tn_wav_(rawTraces, filteredTraces, spikeSites_, spikeTimes, P); fprintf('.');

    if nFet_use >= 2
        spikeSecondarySites = find_site_spk23_(spikeWaveforms, spikeSites_, P);
        spikeWaveforms2 = mn2tn_wav_spk2_(filteredTraces, spikeSecondarySites, spikeTimes, P);
    else
        [spikeSecondarySites, spikeWaveforms2] = deal([]);
    end

    %-----
    % Cancel overlap
    if getOr(P, 'fCancel_overlap', 0)
        try
            [spikeWaveforms, spikeWaveforms2] = cancel_overlap_spk_(spikeWaveforms, spikeWaveforms2, spikeTimes, spikeSites, spikeSecondarySites, siteThresholds, P);
        catch
            fprintf(2, 'fCancel_overlap failed\n');
        end
    end

    spikeTraces = gather_(spikeTraces);
    dialogAssert(nSite_use >0, 'nSites_use = maxSite*2+1 - nSites_ref must be greater than 0');
    switch nFet_use
        case 3
            [spikeSecondarySites, viSite3_spk] = find_site_spk23_(spikeWaveforms, spikeSites_, P); fprintf('.');
            mrFet1 = trWav2fet_(spikeWaveforms, P); fprintf('.');
            mrFet2 = trWav2fet_(spikeWaveforms2, P); fprintf('.');
            mrFet3 = trWav2fet_(mn2tn_wav_spk2_(filteredTraces, viSite3_spk, spikeTimes, P), P); fprintf('.');
            spikeFeatures = permute(cat(3, mrFet1, mrFet2, mrFet3), [1,3,2]); %nSite x nFet x nSpk
            spikePrimarySecondarySites = [spikeSites_(:), spikeSecondarySites(:), viSite3_spk(:)]; %nSpk x nFet
        case 2
            mrFet1 = trWav2fet_(spikeWaveforms, P); fprintf('.');
            mrFet2 = trWav2fet_(spikeWaveforms2, P); fprintf('.');
            spikeFeatures = permute(cat(3, mrFet1, mrFet2), [1,3,2]); %nSite x nFet x nSpk
            spikePrimarySecondarySites = [spikeSites_(:), spikeSecondarySites(:)]; %nSpk x nFet
        case 1
            mrFet1 = trWav2fet_(spikeWaveforms, P); fprintf('.');
            spikeFeatures = permute(mrFet1, [1,3,2]); %nSite x nFet x nSpk
            spikePrimarySecondarySites = [spikeSites_(:)];
        otherwise
            error('wav2spk_: nFet_use must be 1, 2 or 3');
    end

    if nPad_pre > 0
        spikeTimes = spikeTimes - nPad_pre;
    end
    [spikeTimes, spikeFeatures, spikePrimarySecondarySites, spikeWaveforms] = ...
        gather_(spikeTimes, spikeFeatures, spikePrimarySecondarySites, spikeWaveforms);
    useGPU = P.useGPU;

    fprintf('\ttook %0.1fs\n', toc(t_fet));
end % function
