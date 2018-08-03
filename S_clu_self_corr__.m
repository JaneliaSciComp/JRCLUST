%--------------------------------------------------------------------------
function selfcorr = S_clu_self_corr__(S_clu, spikeWaveforms, iClu1, spikeSites)
    % cluster self-correlation. low means bad. return 1-corr score
    MAX_SAMPLE = 4000;
    if nargin<4, spikeSites = get0_('spikeSites'); end

    [viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(S_clu, iClu1, spikeSites);

    viSpk_clu1 = randomSelect_(viSpk_clu1, MAX_SAMPLE);
    % trWav1 = meanSubt_(single(spikeWaveforms(:,:,viSpk_clu1)));
    trWav1 = spikeWaveforms(:,:,viSpk_clu1);
    vrVpp = squeeze_(squeeze_(max(trWav1(:,1,:)) - min(trWav1(:,1,:))));
    % vrVpp = sum(squeeze_(max(spikeWaveforms) - min(trWav1)));
    [~, viSrt] = sort(vrVpp);
    imid = round(numel(viSrt)/2);
    mrWavA = meanSubt_(mean(trWav1(:, :, viSrt(1:imid)), 3));
    mrWavB = meanSubt_(mean(trWav1(:, :, viSrt(imid+1:end)), 3));
    % selfcorr = calc_corr_(mrWavA(:), mrWavB(:));
    % selfcorr = mean(mean(zscore_(mrWavA) .* zscore_(mrWavB)));
    % selfcorr = mean(zscore_(mrWavA(:)) .* zscore_(mrWavB(:)));
    selfcorr = corr_(mrWavA(:), mrWavB(:));
end %func
