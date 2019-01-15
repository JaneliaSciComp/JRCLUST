function mrWavCor = waveformSim(hClust, hCfg, updateMe_)
    %WAVEFORMSIM Compute waveform-based similarity scores for all clusters
    if nargin < 3 || isempty(hClust.mrWavCor)
        updateMe_ = [];
    end

    hCfg.useGPU = 0;

    % not in default param set, but can be overridden if you really want
    fUsePeak2 = hCfg.getOr('fUsePeak2', 0);
    fWaveform_raw = hCfg.getOr('fWavRaw_merge', 1); % revert TW: get_set_(P, 'fWavRaw_merge', 0);
    fZeroStart_raw = hCfg.getOr('fZeroStart_raw', 0);
    fRankCorr_merge = hCfg.getOr('fRankCorr_merge', 0);
    fMode_cor = hCfg.getOr('fMode_cor', 1); % 0: pearson, 1: no mean subt pearson

    if hCfg.verbose
        fprintf('Computing waveform correlation...');
        t1 = tic;
    end

    if fWaveform_raw
        meanWfGlobal = trimRawWaveforms(hClust.meanWfGlobalRaw, hCfg);
    else
        meanWfGlobal = hClust.meanWfGlobal;
    end

    if hCfg.fDrift_merge && fWaveform_raw % only works on raw
        meanWfRawLow = trimRawWaveforms(hClust.meanWfRawLow, hCfg);
        meanWfRawHigh = trimRawWaveforms(hClust.meanWfRawHigh, hCfg);
        meanWfSet = {meanWfGlobal, meanWfRawLow, meanWfRawHigh};
    else
        meanWfSet = {meanWfGlobal};
    end

    if fZeroStart_raw
        meanWfSet = cellfun(@(x) zeroStart(x), meanWfSet, 'UniformOutput', 0);
    end

    if fRankCorr_merge
        meanWfSet = cellfun(@(x) rankorderMatrix(x), meanWfSet, 'UniformOutput', 0);
    end

    if hCfg.nInterp_merge > 1
        meanWfSet = cellfun(@(x) jrclust.utils.interpWindows(x, hCfg.nInterp_merge), meanWfSet, 'UniformOutput', 0);
    end

    if fUsePeak2
        [peaks1, peaks2, peaks3] = hClust.getSecondaryPeaks();
        clusterSites = {peaks1, peaks2, peaks3};
    else
        clusterSites = {hClust.clusterSites};
    end

    % mrWavCor = jrclust.utils.tryGpuArray(zeros(hClust.nClusters), hCfg.useGPU);
    mrWavCor = zeros(hClust.nClusters);

    % shift waveforms up to some multiple (nInterp_merge) of the refractory period
    nShifts = ceil(hCfg.nInterp_merge * hCfg.refracIntms*hCfg.sampleRate/1000);
    % [cviShift1, cviShift2] = jrclust.utils.shiftRange(jrclust.utils.tryGpuArray(int32(size(meanWfGlobal, 1)), hCfg.useGPU), nShifts);
    [cviShift1, cviShift2] = jrclust.utils.shiftRange(int32(size(meanWfGlobal, 1)), nShifts);

    if isempty(updateMe_)
        updateMe = true(hClust.nClusters, 1);
        mrWavCorOld = [];
        fParfor = 1;
    else
        fParfor = 0;
        updateMe = false(hClust.nClusters, 1);
        updateMe(updateMe_) = 1;
        mrWavCorOld = hClust.mrWavCor;
        updateMe((1:hClust.nClusters) > size(mrWavCorOld, 1)) = 1;
    end

    corrData = struct('updateMe', updateMe, ...
                      'cviShift1', {cviShift1}, ...
                      'cviShift2', {cviShift2}, ...
                      'mrWavCorOld', mrWavCorOld, ...
                      'fMode_cor', fMode_cor);

    if fParfor
        try
            parfor iCluster = 1:hClust.nClusters
                pwCor = clu_wavcor_(meanWfSet, clusterSites, hCfg, corrData, iCluster);

                if ~isempty(pwCor)
                    mrWavCor(:, iCluster) = pwCor;
                end
            end
        catch
            fprintf('S_clu_wavcor_: parfor failed. retrying for loop\n');
            fParfor = 0;
        end
    end

    if ~fParfor
        for iCluster = 1:hClust.nClusters
            pwCor = clu_wavcor_(meanWfSet, clusterSites, hCfg, corrData, iCluster);

            if ~isempty(pwCor)
                mrWavCor(:, iCluster) = pwCor;
            end
        end
    end

    mrWavCor = max(mrWavCor, mrWavCor'); % make it symmetric
    mrWavCor(mrWavCor == 0) = nan;
    % mrWavCor = jrclust.utils.tryGather(mrWavCor);

    if hCfg.verbose
        fprintf('\ttook %0.1fs\n', toc(t1));
    end
end