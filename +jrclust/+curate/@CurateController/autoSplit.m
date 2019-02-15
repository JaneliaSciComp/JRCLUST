function autoSplit(obj, multisite)
    %AUTOSPLIT
    if numel(obj.selected) > 1
        return;
    end
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    if obj.hClust.unitCount(obj.selected) < 2
        msgbox('At least two spikes required for splitting');
        return;
    end

    hBox = msgbox('Splitting... (this closes automatically)');
    iCluster = obj.selected(1);
    iSite = obj.hClust.clusterSites(iCluster);

    if multisite
        spikeSites = obj.hCfg.siteNeighbors(1:end - obj.hCfg.nSitesExcl, iSite);
    else
        spikeSites = iSite;
    end

    iSpikes = obj.hClust.spikesByCluster{iCluster};

    sampledSpikes = obj.hClust.getSpikeWindows(iSpikes, spikeSites, 0, 1);
    sampledSpikes = reshape(sampledSpikes, [], size(sampledSpikes, 3));

    % get Vpp of cluster spikes on current site (in FigTime)
    localSpikes = squeeze(obj.hClust.getSpikeWindows(iSpikes, obj.currentSite, 0, 1));
    localVpp = max(localSpikes) - min(localSpikes); % TW calculate amplitudes on the fly

    clusterTimes = obj.hClust.spikeTimes(iSpikes);
    [retained, splitFeatures] = doAutoSplit(sampledSpikes, [clusterTimes localVpp'], 2, obj.hCfg); %TW

    % create split figure
    hFigSplit = jrclust.views.Figure('FigSplit', [.5 0 .5 1], 'Split', 0, 0);
    hFigSplit.addSubplot('pcPlots', 2, 2);

    jrclust.utils.tryClose(hBox);

    while 1
        splitOff = ~retained;
        % plot PC2 vs. PC1
        hFigSplit.subplotApply('pcPlots', 1, @plot, ...
                               splitFeatures(retained, 1), splitFeatures(retained, 2), 'b.', ...
                               splitFeatures(splitOff, 1), splitFeatures(splitOff, 2), 'r.');
        hFigSplit.subplotApply('pcPlots', 1, @xlabel, 'PC 1');
        hFigSplit.subplotApply('pcPlots', 1, @ylabel, 'PC 2');

        % plot PC2 vs. PC3
        hFigSplit.subplotApply('pcPlots', 2, @plot, ...
                               splitFeatures(retained, 3), splitFeatures(retained, 2), 'b.', ...
                               splitFeatures(splitOff, 3), splitFeatures(splitOff, 2), 'r.');
        hFigSplit.subplotApply('pcPlots', 2, @xlabel, 'PC 3');
        hFigSplit.subplotApply('pcPlots', 2, @ylabel, 'PC 2');

        % plot PC3 vs. PC1
        hFigSplit.subplotApply('pcPlots', 3, @plot, ...
                               splitFeatures(retained, 1), splitFeatures(retained, 3), 'b.', ...
                               splitFeatures(splitOff, 1), splitFeatures(splitOff, 3), 'r.');
        hFigSplit.subplotApply('pcPlots', 3, @xlabel, 'PC 1');
        hFigSplit.subplotApply('pcPlots', 3, @ylabel, 'PC 3');

        % plot mean waveforms
        yMin = min(reshape(sampledSpikes, 1, []));
        yMax = max(reshape(sampledSpikes, 1, []));
        meanSamp = jrclust.utils.subsample(1:size(sampledSpikes, 1), 1000);

        hFigSplit.subplotApply('pcPlots', 4, @plot, ...
                               mean(sampledSpikes(meanSamp, retained), 2), 'b');
        hFigSplit.subplotApply('pcPlots', 4, @hold, 'on');
        hFigSplit.subplotApply('pcPlots', 4, @plot, ...
                               mean(sampledSpikes(meanSamp, splitOff), 2), 'r');
        hFigSplit.subplotApply('pcPlots', 4, @ylim, [yMin yMax]);
        hFigSplit.subplotApply('pcPlots', 4, @hold, 'off');
        hFigSplit.subplotApply('pcPlots', 4, @title, 'Mean spike waveforms');

        % ask if we want to split
        dlgAns = questdlg('Split?', 'Confirm split', ...
                          'Yes', 'No', 'Manual', ... % options
                          'No');                     % default

        if strcmp(dlgAns, 'Yes')
            hFigSplit.close();
            break;
        elseif strcmp(dlgAns, 'No')
            hFigSplit.close();
            return;
        else % Manual
            dlgAns = questdlg('Select projection', '', ...
                              'PC2 vs. PC1', 'PC2 vs. PC3', 'PC3 vs. PC1', ... % options
                              'PC2 vs. PC1');                                  % default

            if strcmp(dlgAns,'PC2 vs. PC1')
                spIndex = 1;
                pcX = 1;
                pcY = 2;
            elseif strcmp(dlgAns,'PC2 vs. PC3')
                spIndex = 2;
                pcX = 3;
                pcY = 2;
            else % PC3 vs. PC1
                spIndex = 3;
                pcX = 1;
                pcY = 3;
            end
            % clear axes
            hFigSplit.subplotApply('pcPlots', spIndex, @cla);

            [featuresX, featuresY] = deal(splitFeatures(:, pcX), splitFeatures(:, pcY));
            hFigSplit.subplotApply('pcPlots', spIndex, @plot, featuresX, featuresY, 'k.');

            % user draws a polygon around features to keep
            hPoly = hFigSplit.subplotApply('pcPlots', spIndex, @impoly);
            polyPos = getPosition(hPoly);
            retained = inpolygon(featuresX, featuresY, polyPos(:, 1), polyPos(:, 2));
        end
    end
    hFigSplit.close();

    obj.isWorking = 0;
    obj.splitCluster(iCluster, retained);
end

%% LOCAL FUNCTIONS
function [retained, spikeFeatures] = doAutoSplit(sampledSpikes, spikeFeatures, nSplits, hCfg)
    %DOAUTOSPLIT
    % TODO: ask users number of clusters and split multi-way
    %Make automatic split of clusters using PCA + kmeans clustering
    [~, pcaFeatures, ~] = pca(double(sampledSpikes'), 'Centered', 1, 'NumComponents', 3);
    spikeFeatures = double([spikeFeatures pcaFeatures]);

    if hCfg.getOr('fUseMikeSplit', 0)
        inClust = MikeSplit(sampledSpikes, spikeFeatures, nSplits);
    end

    spikeFeatures = pcaFeatures;
    nSpikes = size(sampledSpikes, 2);

    % ask how many clusters there are
    try
        % kmean clustering into 2
        kmAssigns = kmeans(spikeFeatures, nSplits);
        d12 = madDist(spikeFeatures(kmAssigns == 1, :)', spikeFeatures(kmAssigns == 2, :)');

        if hCfg.verbose
            fprintf('MAD distance: %f\n', d12);
        end

        if hCfg.getOr('fUseMikeSplit', 0)
            retained = inClust;
        else
            retained = (kmAssigns == 1);
        end
    catch % not enough features for k-means to automatically split
        retained = false(nSpikes, 1);
        retained(1:end/2) = 1;
    end
end

function d12 = madDist(features1, features2)
    % distance between two clusters
    if ~ismatrix(features1)
        features1 = reshape(features1, [], size(features1, 3));
    end
    if ~ismatrix(features2)
        features2 = reshape(features2, [], size(features2, 3));
    end

    med1 = median(features1, 2);
    med2 = median(features2, 2);

    medDiffs = med1 - med2;
    norm12 = norm(medDiffs, 2);
    vrFet12_med1 = medDiffs/norm12;

    mad1 = median(abs(vrFet12_med1' * bsxfun(@minus, features1, med1)));
    mad2 = median(abs(vrFet12_med1' * bsxfun(@minus, features2, med2)));

    d12 = norm12 / norm([mad1, mad2], 2);
end