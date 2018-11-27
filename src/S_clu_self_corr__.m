%--------------------------------------------------------------------------
function selfcorr = S_clu_self_corr__(S_clu, tnWav_spk, iClu1, viSite_spk)
    % cluster self-correlation. low means bad. return 1-corr score
    MAX_SAMPLE = 4000;
    if nargin<4, viSite_spk = get0_('viSite_spk'); end

    [viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(S_clu, iClu1, viSite_spk);

    viSpk_clu1 = randomSelect_(viSpk_clu1, MAX_SAMPLE);
    % trWav1 = jrclust.utils.meanSubtract(single(tnWav_spk(:,:,viSpk_clu1)));
    trWav1 = tnWav_spk(:,:,viSpk_clu1);
    vrVpp = squeeze_(squeeze_(max(trWav1(:,1,:)) - min(trWav1(:,1,:))));
    % vrVpp = sum(squeeze_(max(tnWav_spk) - min(trWav1)));
    [~, viSrt] = sort(vrVpp);
    imid = round(numel(viSrt)/2);
    mrWavA = jrclust.utils.meanSubtract(mean(trWav1(:, :, viSrt(1:imid)), 3));
    mrWavB = jrclust.utils.meanSubtract(mean(trWav1(:, :, viSrt(imid+1:end)), 3));
    % selfcorr = calc_corr_(mrWavA(:), mrWavB(:));
    % selfcorr = mean(mean(jrclust.utils.zscore(mrWavA) .* jrclust.utils.zscore(mrWavB)));
    % selfcorr = mean(jrclust.utils.zscore(mrWavA(:)) .* jrclust.utils.zscore(mrWavB(:)));
    selfcorr = corr_(mrWavA(:), mrWavB(:));
end %func
