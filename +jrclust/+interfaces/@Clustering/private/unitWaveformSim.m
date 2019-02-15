function unitSims = unitWaveformSim(meanWfSet, clusterSites, hCfg, scoreData, iCluster)
    %UNITWAVEFORMSIM Mean waveform similarity scores for one cluster against all others
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
        unitSims = [];

        return;
    end

    iNeighbors = hCfg.siteNeighbors(:, iSite);
    
    if fUsePeak2
        compClusters = find(sites == iSite | sites2 == iSite | sites3 == iSite | sites == sites2(iCluster) | sites == sites3(iCluster));
    else
        compClusters = find(ismember(sites, jrclust.utils.findNearbySites(hCfg.siteLoc, iSite, hCfg.evtMergeRad)));
    end

    unitSims = zeros(nClusters, 1, 'single');
    compClusters(compClusters >= iCluster) = []; % symmetric matrix comparison

    if isempty(compClusters)
        return;
    end

    iWaveforms = cellfun(@(x) x(:, iNeighbors, iCluster), meanWfSet, 'UniformOutput', 0);
    jWaveforms = cellfun(@(x) x(:, iNeighbors, compClusters), meanWfSet, 'UniformOutput', 0);

    for j = 1:numel(compClusters)
        jCluster = compClusters(j);

        if ~updateMe(jCluster) && ~updateMe(iCluster)
            unitSims(jCluster) = simScoreOld(jCluster, iCluster);
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
                unitSims(jCluster) = maxCorPearson(iWaveforms_, jWaveforms_, cviShift1, cviShift2, fMode_cor);
            else % dist
                unitSims(jCluster) = minNormDist(iWaveforms_, jWaveforms_);
            end
        end
    end
end

%% LOCAL FUNCTIONS
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