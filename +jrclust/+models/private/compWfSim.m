function simScore = compWfSim(hClust, updateMe_)
    %COMPWFSIM Compute waveform-based similarity scores for all clusters
    % these guys are not in the default param set, but can be overridden if you really want
    fUsePeak2 = hClust.hCfg.getOr('fUsePeak2', false);
    fWaveform_raw = hClust.hCfg.getOr('fWavRaw_merge', true); % revert TW: get_set_(P, 'fWavRaw_merge', 0);
    fZeroStart_raw = hClust.hCfg.getOr('fZeroStart_raw', false);
    fRankCorr_merge = hClust.hCfg.getOr('fRankCorr_merge', false);
    fMode_cor = hClust.hCfg.getOr('fMode_cor', 1); % 0: pearson, 1: no mean subt pearson

    if hClust.hCfg.verbose
        fprintf('Computing waveform correlation...');
        t1 = tic;
    end

    if fWaveform_raw
        meanWfGlobal_ = trimRawWaveforms(hClust.meanWfGlobalRaw, hClust.hCfg);
    else
        meanWfGlobal_ = hClust.meanWfGlobal;
    end

    if hClust.hCfg.fDrift_merge && fWaveform_raw % only works on raw
        meanWfRawLow_ = trimRawWaveforms(hClust.meanWfRawLow, hClust.hCfg);
        meanWfRawHigh_ = trimRawWaveforms(hClust.meanWfRawHigh, hClust.hCfg);
        meanWfSet = {meanWfGlobal_, meanWfRawLow_, meanWfRawHigh_};
    else
        meanWfSet = {meanWfGlobal_};
    end

    if fZeroStart_raw
        meanWfSet = cellfun(@(x) zeroStart(x), meanWfSet, 'UniformOutput', 0);
    end

    if fRankCorr_merge
        meanWfSet = cellfun(@(x) rankorderMatrix(x), meanWfSet, 'UniformOutput', 0);
    end

    if hClust.hCfg.nInterp_merge > 1
        meanWfSet = cellfun(@(x) interpft_(x, hClust.hCfg.nInterp_merge), meanWfSet, 'UniformOutput', 0);
    end

    if fUsePeak2
        [peaks1, peaks2, peaks3] = hClust.getSecondaryPeaks();
        clusterSites_ = {peaks1, peaks2, peaks3};
    else
        clusterSites_ = {hClust.clusterSites};
    end

    % simScore = jrclust.utils.tryGpuArray(zeros(hClust.nClusters), obj.hCfg.useGPU);
    simScore = zeros(hClust.nClusters);

    % shift waveforms up to some multiple (nInterp_merge) of the refractory period
    nShifts = ceil(hClust.hCfg.nInterp_merge * hClust.hCfg.refracIntms*hClust.hCfg.sampleRate/1000);
    % [cviShift1, cviShift2] = jrclust.utils.shiftRange(jrclust.utils.tryGpuArray(int32(size(meanWfGlobal, 1)), obj.hCfg.useGPU), nShifts);
    [cviShift1, cviShift2] = jrclust.utils.shiftRange(int32(size(meanWfGlobal_, 1)), nShifts);

    scoreData = struct('cviShift1', {cviShift1}, ...
                      'cviShift2', {cviShift2}, ...
                      'fMode_cor', fMode_cor);
    if isempty(updateMe_)
        fParfor = true;
        scoreData.updateMe = true(hClust.nClusters, 1);
        scoreData.simScoreOld = [];
    else
        fParfor = false;
        scoreData.updateMe = false(hClust.nClusters, 1);
        scoreData.updateMe(updateMe_) = true;
        scoreData.simScoreOld = hClust.simScore;
        scoreData.updateMe((1:hClust.nClusters) > size(simScoreOld, 1)) = true;
    end

    if fParfor
        try
            parfor iCluster = 1:hClust.nClusters
                pwCor = clusterWaveformSim(meanWfSet, clusterSites_, hClust.hCfg, scoreData, iCluster);

                if ~isempty(pwCor)
                    simScore(:, iCluster) = pwCor;
                end
            end
        catch
            fprintf('computeWaveformSim: parfor failed. retrying for loop\n');
            fParfor = false;
        end
    end

    if ~fParfor
        for iCluster = 1:hClust.nClusters
            pwCor = clusterWaveformSim(meanWfSet, clusterSites_, hClust.hCfg, scoreData, iCluster);

            if ~isempty(pwCor)
                simScore(:, iCluster) = pwCor;
            end
        end
    end

    simScore = max(simScore, simScore'); % make it symmetric
    simScore(simScore == 0) = nan;

    if hClust.hCfg.verbose
        fprintf('\ttook %0.1fs\n', toc(t1));
    end
end

    