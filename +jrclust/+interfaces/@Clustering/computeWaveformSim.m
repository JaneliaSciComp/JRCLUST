function computeWaveformSim(obj, updateMe)
    %COMPUTWAVEFORMSIM Compute waveform-based similarity scores for all clusters
    if nargin < 2 || isempty(obj.waveformSim)
        updateMe = [];
    end

    if isempty(obj.meanWfGlobal)
        return;
    end

    % these guys are not in the default param set, but can be overridden if you really want
    usePeak2 = obj.hCfg.getOr('fUsePeak2', 0);
    useRaw = obj.hCfg.getOr('fWavRaw_merge', 0); % TW
    fZeroStart_raw = obj.hCfg.getOr('fZeroStart_raw', 0);
    fRankCorr_merge = obj.hCfg.getOr('fRankCorr_merge', 0);
    fMode_cor = obj.hCfg.getOr('fMode_cor', 1); % 0: Pearson, 1: non mean-subtracted Pearson

    obj.hCfg.updateLog('wfCorr', 'Computing waveform correlation', 1, 0);

    if useRaw
        meanWfGlobal_ = trimRawWaveforms(obj.meanWfGlobalRaw, obj.hCfg);
    else
        meanWfGlobal_ = obj.meanWfGlobal;
    end

    if obj.hCfg.driftMerge && useRaw % only works on raw
        meanWfRawLow_ = trimRawWaveforms(obj.meanWfRawLow, obj.hCfg);
        meanWfRawHigh_ = trimRawWaveforms(obj.meanWfRawHigh, obj.hCfg);
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

    if obj.hCfg.meanInterpFactor > 1
        meanWfSet = cellfun(@(x) jrclust.utils.interpWindows(x, obj.hCfg.meanInterpFactor), meanWfSet, 'UniformOutput', 0);
    end

    if usePeak2
        [peaks1, peaks2, peaks3] = obj.getSecondaryPeaks();
        clusterSites_ = {peaks1, peaks2, peaks3};
    else
        clusterSites_ = {obj.clusterSites};
    end

    waveformSim = zeros(size(meanWfGlobal_, 3), 'like', obj.waveformSim);

    % shift waveforms up to some multiple (meanInterpFactor) of the refractory period
    nShifts = ceil(obj.hCfg.meanInterpFactor*obj.hCfg.refracInt*obj.hCfg.sampleRate/1000);
    [cviShift1, cviShift2] = jrclust.utils.shiftRange(int32(size(meanWfGlobal_, 1)), nShifts);

    scoreData = struct('cviShift1', {cviShift1}, ...
                      'cviShift2', {cviShift2}, ...
                      'fMode_cor', fMode_cor);

    if isempty(updateMe) % update everything
        useParfor = 1;
        scoreData.updateMe = true(obj.nClusters, 1);
        scoreData.simScoreOld = [];
    else
        useParfor = 0;
        scoreData.updateMe = false(obj.nClusters, 1);
        scoreData.updateMe(updateMe) = 1;
        scoreData.simScoreOld = obj.waveformSim;
        scoreData.updateMe((1:obj.nClusters) > size(scoreData.simScoreOld, 1)) = 1;
    end

    if useParfor
        % avoid sending the entire hCfg object out to workers
        cfgSub = struct('siteNeighbors', obj.hCfg.siteNeighbors, ...
                        'siteLoc', obj.hCfg.siteLoc, ...
                        'evtMergeRad', obj.hCfg.evtMergeRad, ...
                        'autoMergeBy', obj.hCfg.autoMergeBy);
        nClusters_ = obj.nClusters;
        try
            parfor iCluster = 1:nClusters_
                pwCor = unitWaveformSim(meanWfSet, clusterSites_, cfgSub, scoreData, iCluster);

                if ~isempty(pwCor)
                    waveformSim(:, iCluster) = pwCor;
                end
            end
        catch ME
            warning('computeWaveformSim: parfor failed');
            useParfor = 0;
        end
    end

    if ~useParfor
        for iCluster = 1:obj.nClusters
            pwCor = unitWaveformSim(meanWfSet, clusterSites_, obj.hCfg, scoreData, iCluster);

            if ~isempty(pwCor)
                waveformSim(:, iCluster) = pwCor;
            end
        end
    end

    waveformSim = max(waveformSim, waveformSim'); % make it symmetric

    obj.hCfg.updateLog('wfCorr', 'Finished computing waveform correlation', 0, 1);

    if isempty(updateMe)
        waveformSim = jrclust.utils.setDiag(waveformSim, obj.computeSelfSim());
    else
        % carry over the old diagonal
        waveformSim = jrclust.utils.setDiag(waveformSim, diag(obj.waveformSim));

        for ii = 1:numel(updateMe)
            iCluster = updateMe(ii);
            waveformSim(iCluster, iCluster) = obj.computeSelfSim(iCluster);
        end
    end

    %waveformSim(waveformSim == 0) = nan;
    obj.waveformSim = waveformSim;
end

%% LOCAL FUNCTIONS
function mi = rankorderMatrix(meanWf)
    %RANKORDERMATRIX 
    if isrow(meanWf)
        meanWf = meanWf';
    end

    shape = size(meanWf);
    if numel(shape) == 3
        meanWf = meanWf(:); % global order
    end

    mi = zeros(size(meanWf));
    for iCol = 1:size(meanWf, 2)
        col = meanWf(:, iCol);
        pos = find(col > 0);

        if ~isempty(pos)
            [~, argsort] = sort(col(pos), 'ascend');
            mi(pos(argsort), iCol) = 1:numel(pos);
        end

        neg = find(col < 0);
        if ~isempty(neg)
            [~, argsort] = sort(col(neg), 'descend');
            mi(neg(argsort), iCol) = -(1:numel(neg));
        end
    end

    if numel(shape) == 3
        mi = reshape(mi, shape);
    end
end

function meanWf = trimRawWaveforms(meanWf, hCfg)
    %TRIMRAWWAVEFORMS Pare raw waveforms down to the size of filtered waveforms
    nSamplesRaw = diff(hCfg.evtWindowRawSamp) + 1;
    spkLimMerge = round(hCfg.evtWindowSamp * hCfg.evtWindowMergeFactor);
    nSamplesRawMerge = diff(spkLimMerge) + 1;

    if size(meanWf, 1) <= nSamplesRawMerge
        return;
    end

    limits = [spkLimMerge(1) - hCfg.evtWindowRawSamp(1) + 1, ...
              nSamplesRaw - (hCfg.evtWindowRawSamp(2) - spkLimMerge(2))];

    meanWf = meanWf(limits(1):limits(2), :, :);
    meanWf = jrclust.utils.meanSubtract(meanWf);
end

function meanWf = zeroStart(meanWf)
    %ZEROSTART Subtract the first row from all rows
    shape = size(meanWf);
    if numel(shape) ~= 2
        meanWf = reshape(meanWf, shape(1), []);
    end

    meanWf = bsxfun(@minus, meanWf, meanWf(1, :));
    if numel(shape) ~= 2
        meanWf = reshape(meanWf, shape);
    end
end