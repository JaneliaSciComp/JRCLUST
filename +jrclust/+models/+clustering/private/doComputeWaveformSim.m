function simScore = doComputeWaveformSim(hClust, updateMe)
    %COMPWFSIM Compute waveform-based similarity scores for all clusters
    % these guys are not in the default param set, but can be overridden if you really want
    fUsePeak2 = hClust.hCfg.getOr('fUsePeak2', 0);
    fWaveform_raw = hClust.hCfg.getOr('fWavRaw_merge', 0); % TW
    fZeroStart_raw = hClust.hCfg.getOr('fZeroStart_raw', 0);
    fRankCorr_merge = hClust.hCfg.getOr('fRankCorr_merge', 0);
    fMode_cor = hClust.hCfg.getOr('fMode_cor', 1); % 0: Pearson, 1: non mean-subtracted Pearson

    if hClust.hCfg.verbose
        fprintf('Computing waveform correlation...');
        t1 = tic;
    end

    if fWaveform_raw
        meanWfGlobal_ = trimRawWaveforms(hClust.meanWfGlobalRaw, hClust.hCfg);
    else
        meanWfGlobal_ = hClust.meanWfGlobal;
    end

    if hClust.hCfg.driftMerge && fWaveform_raw % only works on raw
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

    if hClust.hCfg.meanInterpFactor > 1
        meanWfSet = cellfun(@(x) jrclust.utils.interpWindows(x, hClust.hCfg.meanInterpFactor), meanWfSet, 'UniformOutput', 0);
    end

    if fUsePeak2
        [peaks1, peaks2, peaks3] = hClust.getSecondaryPeaks();
        clusterSites_ = {peaks1, peaks2, peaks3};
    else
        clusterSites_ = {hClust.clusterSites};
    end

    simScore = zeros(size(meanWfGlobal_, 3), 'like', hClust.simScore);

    % shift waveforms up to some multiple (meanInterpFactor) of the refractory period
    nShifts = ceil(hClust.hCfg.meanInterpFactor*hClust.hCfg.refracInt*hClust.hCfg.sampleRate/1000);
    [cviShift1, cviShift2] = jrclust.utils.shiftRange(int32(size(meanWfGlobal_, 1)), nShifts);

    scoreData = struct('cviShift1', {cviShift1}, ...
                      'cviShift2', {cviShift2}, ...
                      'fMode_cor', fMode_cor);
    if isempty(updateMe)
        useParfor = 0;
        scoreData.updateMe = true(hClust.nClusters, 1);
        scoreData.simScoreOld = [];
    else
        useParfor = 0;
        scoreData.updateMe = false(hClust.nClusters, 1);
        scoreData.updateMe(updateMe) = 1;
        scoreData.simScoreOld = hClust.simScore;
        scoreData.updateMe((1:hClust.nClusters) > size(scoreData.simScoreOld, 1)) = 1;
    end

    if useParfor
        % avoid sending the entire hCfg object out to workers
        cfgSub = struct('siteNeighbors', hClust.hCfg.siteNeighbors, ...
                        'siteLoc', hClust.hCfg.siteLoc, ...
                        'evtDetectRad', hClust.hCfg.evtDetectRad, ...
                        'autoMergeBy', hClust.hCfg.autoMergeBy);
        try
            parfor iCluster = 1:hClust.nClusters
                pwCor = clusterWaveformSim(meanWfSet, clusterSites_, cfgSub, scoreData, iCluster);

                if ~isempty(pwCor)
                    simScore(:, iCluster) = pwCor;
                end
            end
        catch ME
            warning('computeWaveformSim: parfor failed: %s', ME.messsage);
            useParfor = 0;
        end
    end

    if ~useParfor
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

%% LOCAL FUNCTIONS
function corScores = clusterWaveformSim(meanWfSet, clusterSites, hCfg, scoreData, iCluster)
    %CLUSTERWAVEFORMSIM Mean waveform similarity scores for one cluster against all others
    updateMe = scoreData.updateMe;
    cviShift1 = scoreData.cviShift1;
    cviShift2 = scoreData.cviShift2;
    simScoreOld = scoreData.simScoreOld;
    fMode_cor = scoreData.fMode_cor;

    if numel(clusterSites) == 1
        sites = clusterSites{1};
        fUsePeak2 = 0;
    else
        sites = clusterSites{1};
        sites2 = clusterSites{2};
        sites3 = clusterSites{3};
        fUsePeak2 = 1;
    end

    nClusters = numel(sites);
    iSite = sites(iCluster);

    if iSite == 0 || isnan(iSite)
        corScores = [];

        return;
    end

    iNeighbors = hCfg.siteNeighbors(:, iSite);
    
    if fUsePeak2
        compClusters = find(sites == iSite | sites2 == iSite | sites3 == iSite | sites == sites2(iCluster) | sites == sites3(iCluster));
    else
        compClusters = find(ismember(sites, jrclust.utils.findNearbySites(hCfg.siteLoc, iSite, hCfg.evtMergeRad)));
    end

    corScores = zeros(nClusters, 1, 'single');
    compClusters(compClusters >= iCluster) = []; % symmetric matrix comparison

    if isempty(compClusters)
        return;
    end

    iWaveforms = cellfun(@(x) x(:, iNeighbors, iCluster), meanWfSet, 'UniformOutput', 0);
    jWaveforms = cellfun(@(x) x(:, iNeighbors, compClusters), meanWfSet, 'UniformOutput', 0);

    for j = 1:numel(compClusters)
        jCluster = compClusters(j);

        if ~updateMe(jCluster) && ~updateMe(iCluster)
            corScores(jCluster) = simScoreOld(jCluster, iCluster);
        else
            jSite = sites(jCluster);
            if jSite == 0 || isnan(jSite)
                continue;
            end

            if jSite == iSite % neighbors are the same, can compare all sites
                iWaveforms_ = iWaveforms;
                jWaveforms_ = cellfun(@(x) x(:, :, j), jWaveforms, 'UniformOutput', 0);
            else % find overlap of neighbors and compare these
                jNeighbors = hCfg.siteNeighbors(:, jSite);
                ijOverlap = find(ismember(iNeighbors, jNeighbors));

                if isempty(ijOverlap) % nothing to compare
                    continue;
                end

                iWaveforms_ = cellfun(@(x) x(:, ijOverlap), iWaveforms, 'UniformOutput', 0);
                jWaveforms_ = cellfun(@(x) x(:, ijOverlap, j), jWaveforms, 'UniformOutput', 0);
            end

            if strcmp(hCfg.autoMergeBy, 'pearson')
                corScores(jCluster) = maxCorPearson(iWaveforms_, jWaveforms_, cviShift1, cviShift2, fMode_cor);
            else % dist
                corScores(jCluster) = minNormDist(iWaveforms_, jWaveforms_);
            end
        end
    end
end

function maxCor = maxCorPearson(iWaveforms, jWaveforms, shifts1, shifts2, fMode_cor)
    assert(numel(iWaveforms) == numel(jWaveforms), 'maxCorPearson: numel must be the same');

    if numel(iWaveforms) == 1
        maxCor = max(shiftPearson(iWaveforms{1}, jWaveforms{1}, shifts1, shifts2));
    else
        tr1 = cat(3, iWaveforms{:}); %nT x nC x nDrifts
        tr2 = cat(3, jWaveforms{:});

        nDrifts = numel(iWaveforms);
        nShifts = numel(shifts1);
        corScores = zeros(1, nShifts);

        for iShift = 1:nShifts
            mr1_ = reshape(tr1(shifts1{iShift}, :, :), [], nDrifts);
            mr2_ = reshape(tr2(shifts2{iShift}, :, :), [], nDrifts);

            if fMode_cor == 0 % center the waveforms
                mr1_ = bsxfun(@minus, mr1_, mean(mr1_));
                mr2_ = bsxfun(@minus, mr2_, mean(mr2_));
                mr1_ = bsxfun(@rdivide, mr1_, sqrt(sum(mr1_.^2)));
                mr2_ = bsxfun(@rdivide, mr2_, sqrt(sum(mr2_.^2)));

                corScores(iShift) = max(max(mr1_' * mr2_));
            else
                mr12_ = (mr1_' * mr2_) ./ sqrt(sum(mr1_.^2)' * sum(mr2_.^2));
                corScores(iShift) = max(mr12_(:));
            end
        end

        maxCor = max(corScores);
    end
end

function minDist = minNormDist(iWaveforms, jWaveforms) % TW
    jrclust.utils.dlgAssert(numel(iWaveforms) == numel(jWaveforms), 'minNormDist: numel must be the same');

    if numel(iWaveforms) == 1
        x = reshape(iWaveforms{1}, [], 1);
        y = reshape(jWaveforms{1}, [], 1);

        minDist = 1 - norm(x - y)/max(norm(x), norm(y));
    else
        tr1 = cat(3, iWaveforms{:}); %nT x nC x nDrifts
        tr2 = cat(3, jWaveforms{:});
        dists = zeros(size(tr1, 3), 1);

        for j = 1:size(tr1, 3)
            x = reshape(tr1(:, :, j), [], 1);
            y = reshape(tr2(:, :, j), [], 1);

            dists(j) = 1 - norm(x - y)/max(norm(x), norm(y));
        end

        minDist = max(dists);
    end
end

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

function scores = shiftPearson(X, Y, shift1, shift2)
    % TODO: use shift matrix?  (https://en.wikipedia.org/wiki/Shift_matrix)
    %   maybe just use shifted indices
    scores = zeros(size(shift1), 'like', X);

    for iShift = 1:numel(scores)
        vX = X(shift1{iShift}, :);
        vY = Y(shift2{iShift}, :);

        % Pearson correlation coefficient for shifted X and Y
        r = corrcoef(vX(:), vY(:)); % TODO: get p-values?
        scores(iShift) = r(1, 2);
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