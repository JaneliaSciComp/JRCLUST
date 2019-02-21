function spikeData = extractFeatures(obj, spikeData)
    %EXTRACTFEATURES Extract spike waveforms and build a spike table
    spikesFilt = spikeData.spikesFilt;
    spikesFilt2 = spikeData.spikesFilt2;
    spikesFilt3 = spikeData.spikesFilt3;

    obj.hCfg.updateLog('extractFeatures', 'Extracting features', 1, 0);

    % spikesRaw, spikesFilt are nSamples x nSites x nSpikes
    %[windowsRaw, windowsFilt, spikeTimes] = obj.samplesToWindows(samplesRaw, samplesFilt, spikeTimes, spikeSites_);

    % get secondary/tertiary peaks and extract windows from them to compute features
%     if obj.hCfg.nPeaksFeatures >= 2
%         [spikeSites2, spikeSites3] = obj.findSecondaryPeaks(windowsFilt, spikeSites_);
%         windowsFilt2 = obj.samplesToWindows2(samplesFilt, spikeSites2, spikeTimes);
% 
%         if obj.hCfg.nPeaksFeatures == 3
%             windowsFilt3 = obj.samplesToWindows2(samplesFilt, spikeSites3, spikeTimes);
%         end
%     else
%         [spikeSites2, windowsFilt2] = deal([]);
%     end

    if obj.hCfg.nPeaksFeatures == 1
        features1 = jrclust.features.computeFeatures(spikesFilt, obj.hCfg);
        % if obj.hCfg.verbose
        %     fprintf('.');
        % end

        spikeFeatures = permute(features1, [1, 3, 2]); % nSites x nFeatures x nSpikes
    elseif obj.hCfg.nPeaksFeatures == 2
        features1 = jrclust.features.computeFeatures(spikesFilt, obj.hCfg);
        % if obj.hCfg.verbose
        %     fprintf('.');
        % end

        features2 = jrclust.features.computeFeatures(spikesFilt2, obj.hCfg);
        % if obj.hCfg.verbose
        %     fprintf('.');
        % end

        spikeFeatures = permute(cat(3, features1, features2), [1, 3, 2]); % nSites x nFeatures x nSpikes
    else % obj.hCfg.nPeaksFeatures == 3
        features1 = jrclust.features.computeFeatures(spikesFilt, obj.hCfg);
        % if obj.hCfg.verbose
        %     fprintf('.');
        % end

        features2 = jrclust.features.computeFeatures(spikesFilt2, obj.hCfg);
        % if obj.hCfg.verbose
        %     fprintf('.');
        % end

        features3 = jrclust.features.computeFeatures(spikesFilt3, obj.hCfg);
        % if obj.hCfg.verbose
        %     fprintf('.');
        % end

        spikeFeatures = permute(cat(3, features1, features2, features3), [1, 3, 2]); % nSites x nFeatures x nSpikes
    end

    spikeData.spikeFeatures = jrclust.utils.tryGather(spikeFeatures);
    obj.hCfg.updateLog('extractFeatures', 'Finished extracting features', 0, 1);
end
