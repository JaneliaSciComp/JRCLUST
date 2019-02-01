function [spikeTimes, spikeAmps, spikeSites] = detectPeaks(samplesIn, siteThresh, keepMe, hCfg)
    %DETECTPEAKS Detect peaks for each site, 
    nSites = size(samplesIn, 2);
    [spikesBySite, ampsBySite] = deal(cell(nSites, 1));

    if hCfg.verbose
        fprintf('\tDetecting spikes from each channel.\n\t\t');
        tDetect = tic;
    end

    for iSite = 1:nSites
        % find spikes
        [peakLocs, peaks] = detectPeaksSite(samplesIn(:, iSite), siteThresh(iSite), hCfg);
        if hCfg.verbose
            fprintf('.');
        end

        if isempty(keepMe)
            spikesBySite{iSite} = peakLocs;
            ampsBySite{iSite} = peaks;
        else % reject global mean
            spikesBySite{iSite} = peakLocs(keepMe(peakLocs));
            ampsBySite{iSite}   = peaks(keepMe(peakLocs));
        end
    end

    if hCfg.verbose
        nSpikes = sum(cellfun(@numel, spikesBySite));
        fprintf('\n\t\tDetected %d spikes from %d sites; took %0.1fs.\n', nSpikes, nSites, toc(tDetect));
    end

    % Group spiking events using vrWav_mean1. already sorted by time
    if hCfg.getOr('fMerge_spk', 1)
        if hCfg.verbose
            fprintf('\tMerging spikes...');
            t2 = tic;
        end

        [spikeTimes, spikeAmps, spikeSites] = mergePeaks(spikesBySite, ampsBySite, hCfg);

        if hCfg.verbose
            fprintf('\t%d spiking events found; took %0.1fs\n', numel(spikeSites), toc(t2));
        end
    else
        spikeTimes = jrclust.utils.neCell2mat(spikesBySite);
        spikeAmps = jrclust.utils.neCell2mat(ampsBySite);

        % generate a bunch of 1's, 2's, ..., nSites's to sort later
        nSpikesSite = cellfun(@(ss) numel(ss), spikesBySite);
        siteOnes = cell(numel(spikesBySite), 1);
        for iSite = 1:numel(spikesBySite)
            siteOnes{iSite} = iSite * ones(nSpikesSite(iSite), 1);
        end
        spikeSites = jrclust.utils.neCell2mat(siteOnes);

        % sort by time
        [spikeTimes, argsort] = sort(spikeTimes, 'ascend');
        spikeAmps = spikeAmps(argsort);
        spikeSites = spikeSites(argsort);
    end

    spikeAmps = jrclust.utils.tryGather(spikeAmps);

    % Group all sites in the same shank
    if hCfg.groupShank
        spikeSites = groupByShank(spikeSites, hCfg); % change the site location to the shank center
    end
end

%% LOCAL FUNCTIONS
function [spikeSites] = groupByShank(spikeSites, hCfg)
    %GROUPBYSHANK Group all spike sites by shank
    site2site = zeros([hCfg.nSites, 1], 'like', spikeSites);

    % remap
    [~, ia, ic] = unique(hCfg.shankMap);
    site2site(hCfg.siteMap) = ia(ic);
    site2site = site2site(site2site > 0);
    spikeSites = site2site(spikeSites);
end
