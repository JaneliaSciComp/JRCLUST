%--------------------------------------------------------------------------
function [S_gt, tnWav_spk, tnWav_raw] = gt2spk_(S_gt, P, snr_thresh)
    % convert ground truth to spike waveforms
    fSubtract_nmean = 0;
    % fSubtract_ref = 1;
    % P.useGPU = 0;
    MAX_SAMPLE = 1000; %mean calc
    nSubsample_clu = get_set_(P, 'nSubsample_clu', 8);

    if nargin<3, snr_thresh = []; end

    fProcessRaw = (nargout == 3);
    % fProcessRaw = 0;
    t1 = tic;
    fprintf('Computing ground truth units...\n');
    viClu = int32(S_gt.viClu);
    viTime_spk = int32(S_gt.viTime);
    nSites = numel(P.viSite2Chan);

    % load entire raw waveform to memory
    % P = setfield(P, 'vcFilter', 'bandpass'); % use bandpass for GT evaluation
    if ~P.fTranspose_bin
        [mnWav, vrWav_mean] = load_file_(P.vcFile, [], P);
        mnWav = fft_clean_(mnWav, P);
        if fProcessRaw
            tnWav_raw = permute(mn2tn_gpu_(mnWav, P.spkLim_raw, viTime_spk), [1,3,2]);
        end
        [mnWav, ~] = filt_car_(mnWav, P);
        %     mnWav = mnWav_filt_(mnWav, P); % Apply filtering in RAM
        if fSubtract_nmean % Apply nmean CAR to ground truth spikes (previous standard)
            P1=P; P1.vcCommonRef = 'nmean'; mnWav = wav_car_(mnWav, P1);
        end
        tnWav_spk = permute(mn2tn_gpu_(mnWav, P.spkLim, viTime_spk), [1,3,2]);
        vrVrms_site = gather_(mr2rms_(mnWav, 1e5));
        clear mnWav;
    else % real ground truth: must block load and filter.
        viClu = viClu(1:nSubsample_clu:end);
        viTime_spk = viTime_spk(1:nSubsample_clu:end);
        [tnWav_spk, vrVrms_site] = file2spk_gt_(P, viTime_spk);
        snr_thresh = []; %no SNR threshold
    end

    % trim excess spikes
    % nSpk = size(tnWav_spk,3);
    % [viClu, viTime_spk, S_gt.viTime, S_gt.viClu] = multifun(@(x)x(1:nSpk), viClu, viTime_spk, S_gt.viTime, S_gt.viClu);

    % determine mean spikes
    nClu = max(viClu);
    trWav_clu = zeros(size(tnWav_spk,1), nSites, nClu, 'single');
    if fProcessRaw
        trWav_raw_clu = zeros(size(tnWav_raw,1), nSites, nClu, 'single');
    else
        trWav_raw_clu = [];
    end
    cviSpk_clu = arrayfun(@(iClu)int32(find(viClu == iClu)), 1:nClu, 'UniformOutput', 0);
    for iClu=1:nClu
        viSpk_clu1 = cviSpk_clu{iClu};
        %     viSpk_clu1 = viSpk_clu1(
        viSpk1 = subsample_vr_(viSpk_clu1, MAX_SAMPLE);
        if isempty(viSpk1), continue; end
        try
            trWav_clu(:,:,iClu) = mean(tnWav_spk(:,:,viSpk1), 3); %multiply by scaling factor?
        catch
            ;
        end
        %     trWav_clu(:,:,iClu) = fft_align_mean_(single(tnWav_spk(:,:,viSpk1)));
        if fProcessRaw
            trWav_raw_clu(:,:,iClu) = mean_tnWav_raw_(tnWav_raw(:,:,viSpk1), P);
        end
        fprintf('.');
    end

    % if fSubtract_ref
    % % apply surrounding reference subtraction
    %     trWav_clu = spkwav_car_(trWav_clu, spkwav_car_init_(P));
    % end

    % Find center location and spike SNR
    mrVmin_clu = shiftdim(min(trWav_clu));
    [vrVmin_clu, viSite_clu] = min(mrVmin_clu); %center sites

    % perform CAR after centering at the center site
    mrVmin_clu = squeeze_(min(trWav_clu));
    [vrVmin_clu, ~] = min(mrVmin_clu);

    % cluster specifications
    vrVmin_clu = abs(vrVmin_clu);
    vrSnr_clu = (vrVmin_clu ./ vrVrms_site(viSite_clu))';
    vrThresh_site = vrVrms_site * P.qqFactor;
    vnSite_clu = sum(bsxfun(@lt, mrVmin_clu, -vrThresh_site(viSite_clu)));

    if ~isempty(snr_thresh)
        viClu_keep = find(abs(vrSnr_clu) > snr_thresh);
        [trWav_clu, trWav_raw_clu] = multifun_(@(x)x(:,:,viClu_keep), trWav_clu, trWav_raw_clu);
        [viSite_clu, vrVmin_clu, vrSnr_clu, cviSpk_clu, vnSite_clu, vrVrms_site] = ...
        multifun_(@(x)x(viClu_keep), viSite_clu, vrVmin_clu, vrSnr_clu, cviSpk_clu, vnSite_clu, vrVrms_site);
        vlSpk_keep = ismember(viClu, viClu_keep);
        [S_gt.viClu, S_gt.viTime] = multifun_(@(x)x(vlSpk_keep), S_gt.viClu, S_gt.viTime);
        viMap = 1:max(viClu_keep);
        viMap(viClu_keep) = 1:numel(viClu_keep);
        S_gt.viClu = viMap(S_gt.viClu); % compact, no-gap
    else
        viClu_keep = [];
    end

    miSites_clu = P.miSites(:,viSite_clu);
    S_gt = struct_add_(S_gt, trWav_clu, trWav_raw_clu, viSite_clu, vrVmin_clu, ...
    vrSnr_clu, cviSpk_clu, vnSite_clu, vrVrms_site, miSites_clu, viClu_keep);
    fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func
