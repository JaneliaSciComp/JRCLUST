%function [vrFet1, vrTime1, vcYlabel, viSpk1] = getDispFeaturesSite(iSite, iCluster, S0)
function [dispFeatures, spikeTimesSecs, yLabel, dispSpikes] = getDispFeaturesSite(hClust, iSite, iCluster)
    %GETDISPFEATURESSITE Compute features to display on a given site,
    %optionally for a given cluster
    if nargin < 3 % get features off of background spikes
        iCluster = [];
    end

    hCfg = hClust.hCfg;
    [dispFeatures, dispSpikes] = getDispFeaturesCluster(hClust, iCluster, iSite);
    spikeTimesSecs = double(hClust.spikeTimes(dispSpikes))/hCfg.sampleRate;

    if strcmp(hCfg.dispFeature, 'vpp')
        yLabel = sprintf('Site %d (\\mu Vpp)', iSite);

    else
        yLabel = sprintf('Site %d (%s)', iSite, hCfg.dispFeature);
    end
end
