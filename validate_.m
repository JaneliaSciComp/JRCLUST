%--------------------------------------------------------------------------
function validate_(P)
    persistent S_gt vcFile_prm_ % spikeWaveforms tnWav_gt
    % S0 = load_cached_(P, 0);
    fMergeCheck = 0; %kilosort-style validation
    % snr_thresh_score = 10;
    snr_thresh_stat = getOr(P, 'snr_thresh_gt', 7);

    S0 = load_(strrep(P.paramFile, '.prm', '_jrc.mat'));
    % if isempty(spikeWaveforms)
    %     spikeWaveforms = load_bin_(strrep(P.paramFile, '.prm', '_waveforms.bin'), P.dataType, S0.waveformDims);
    % end
    S_clu = S0.S_clu;

    % Load ground truth file
    if ~exist(P.groundTruthFile, 'file'), P.groundTruthFile = subsFileExt_(P.vcFile, '_gt.mat'); end
    if ~strcmpi(vcFile_prm_, P.paramFile), vcFile_prm_=[]; end
    if isempty(vcFile_prm_) || isempty(S_gt) % cache ground truth result
        vcFile_prm_ = P.paramFile;
        groundTruthFile1 = strrep(P.paramFile, '.prm', '_gt1.mat');
        if fileExists(groundTruthFile1)
            S_gt = load(groundTruthFile1);
        else
            S_gt0 = load_gt_(P.groundTruthFile, P);
            if isempty(S_gt0), fprintf(2, 'Groundtruth does not exist. Run "jrclust import" to create a groundtruth file.\n'); return; end
            %[S_gt, tnWav_gt] = gt2spk_(S_gt, P, snr_thresh_stat);
            S_gt = gt2spk_(S_gt0, P, snr_thresh_stat);
            struct_save_(S_gt, groundTruthFile1);
        end
    end
    S_score = struct(...
    'vrVmin_gt', S_gt.vrVmin_clu, 'vnSite_gt', S_gt.vnSite_clu, ...
    'vrSnr_gt', S_gt.vrSnr_clu, 'vrSnr_min_gt', S_gt.vrSnr_clu, ...
    'trWav_gt', S_gt.trWav_clu, 'viSite_gt', S_gt.clusterSites);
    S_score.cviSpk_gt = S_gt.spikesByCluster;

    % Compare S_clu with S_gt
    nSamples_jitter = round(P.sampleRateHz / 1000); %1 ms jitter
    fprintf('verifying cluster...\n');
    [mrMiss, mrFp, vnCluGt, miCluMatch, S_score_clu] = ...
    clusterVerify(S_gt.viClu, S_gt.viTime, S_clu.spikeClusters, S0.spikeTimes, nSamples_jitter); %S_gt.viTime
    % viClu_spk = S_score_clu.cviHit_gt
    if fMergeCheck, compareClustering2_(S_gt.viClu, S_gt.viTime, S_clu.spikeClusters+1, S0.spikeTimes); end

    Sgt = S_gt; %backward compatibility
    S_score = struct_add_(S_score, mrMiss, mrFp, vnCluGt, miCluMatch, P, Sgt, S_score_clu);
    S_score.cviTime_clu = S_clu.spikesByCluster(S_score_clu.viCluMatch)';
    S_score.vrVrms_site = single(S0.vrThresh_site) / S0.P.qqFactor;

    fprintf('SNR_gt (Vp/Vrms): %s\n', sprintf('%0.1f ', S_score.vrSnr_gt));
    fprintf('nSites>thresh (GT): %s\n', sprintf('%d ', S_score.vnSite_gt));

    write_struct_(strrep(P.paramFile, '.prm', '_score.mat'), S_score);

    setUserData(S_score);
    assignWorkspace_(S_score); %put in workspace

    figure; set(gcf,'Name',P.paramFile);
    subplot 121; plot_cdf_(S_score.S_score_clu.vrFp); hold on; plot_cdf_(S_score.S_score_clu.vrMiss);
    legend({'False Positives', 'False Negatives'}); ylabel('CDF'); grid on; xlabel('Cluster count');

    subplot 122; hold on;
    plot(S_score.vrSnr_min_gt, S_score.S_score_clu.vrFp, 'b.', S_score.vrSnr_min_gt, S_score.S_score_clu.vrMiss, 'r.');
    legend({'False Positives', 'False Negatives'}); ylabel('score'); grid on; xlabel('SNR (Vp/Vrms)');

    vnSpk_gt = cellfun(@numel, S_score.S_score_clu.cviSpk_gt_hit) + cellfun(@numel, S_score.S_score_clu.cviSpk_gt_miss);
    disp_score_(S_score.vrSnr_min_gt, S_score.S_score_clu.vrFp, S_score.S_score_clu.vrMiss, ...
    S_score.S_score_clu.vrAccuracy, S_score.vnSite_gt, vnSpk_gt, 1);
end % function
