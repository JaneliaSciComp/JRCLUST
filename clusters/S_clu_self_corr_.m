%--------------------------------------------------------------------------
function selfcorr = S_clu_self_corr_(S_clu, iClu1, S0)
    % plot top half vs bottom half correlation. sum of vpp
    if nargin<2, iClu1 = []; end
    if nargin<3, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    [spikeSites, P] = deal(S0.spikeSites, S0.P);
    spikeTraces = getSpikeWaveforms(P, getOr(P, 'fWavRaw_merge', 1));

    if isempty(iClu1)
        fprintf('Computing self correlation\n\t'); t1=tic;
        selfcorr = zeros(1, S_clu.nClusters);
        for iClu=1:S_clu.nClusters
            selfcorr(iClu) = S_clu_self_corr__(S_clu, spikeTraces, iClu, spikeSites);
            fprintf('.');
        end
        fprintf('\n\ttook %0.1fs\n', toc(t1));
    else
        selfcorr = S_clu_self_corr__(S_clu, spikeTraces, iClu1, spikeSites);
    end
end % function

%--------------------------------------------------------------------------
function selfcorr = S_clu_self_corr__(S_clu, spikeWaveforms, iClu1, spikeSites)
    % cluster self-correlation. low means bad. return 1-corr score
    MAX_SAMPLE = 4000;
    if nargin<4, spikeSites = get0_('spikeSites'); end

    [viSpk_clu1, viiSpk_clu1] = getClusterCenteredSpikes(S_clu, iClu1, spikeSites);

    viSpk_clu1 = randomSubsample(viSpk_clu1, MAX_SAMPLE);
    % trWav1 = meanSubtract(single(spikeWaveforms(:,:,viSpk_clu1)));
    trWav1 = spikeWaveforms(:,:,viSpk_clu1);
    vrVpp = squeeze_(squeeze_(max(trWav1(:,1,:)) - min(trWav1(:,1,:))));
    % vrVpp = sum(squeeze_(max(spikeWaveforms) - min(trWav1)));
    [~, viSrt] = sort(vrVpp);
    imid = round(numel(viSrt)/2);
    mrWavA = meanSubtract(mean(trWav1(:, :, viSrt(1:imid)), 3));
    mrWavB = meanSubtract(mean(trWav1(:, :, viSrt(imid+1:end)), 3));
    % selfcorr = calc_corr_(mrWavA(:), mrWavB(:));
    % selfcorr = mean(mean(zscore_(mrWavA) .* zscore_(mrWavB)));
    % selfcorr = mean(zscore_(mrWavA(:)) .* zscore_(mrWavB(:)));
    selfcorr = corr_(mrWavA(:), mrWavB(:));
end % function
