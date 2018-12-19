function [mrPc1, mrPc2, mrPv1, mrPv2] = pca_pc_spk_(hClust, sampledSpikes, iSite, mrPv1, mrPv2)
    %GETSPIKEPC Project selected spikes onto first and second principal components
    nSites = numel(iSite);
    sampledWindows = permute(jrclust.utils.getSampledWindows(hClust, sampledSpikes, iSite, false), [1, 3, 2]); % nSamples x nSpikes x nSites
    if nargin < 4
        [mrPv1, mrPv2] = pca_pv_spk_(sampledSpikes, iSite, sampledWindows);
    end

    [nSamples, nSpikes, ~] = size(sampledWindows);
    [mrPc1, mrPc2] = deal(zeros(nSpikes, nSites, 'single'));

    try
        for jSite = 1:nSites
            jWindow = jrclust.utils.meanSubtract(single(sampledWindows(:, :, jSite)));
            mrPc1(:, jSite) = (mrPv1(:, jSite)' * jWindow)';
            mrPc2(:, jSite) = (mrPv2(:, jSite)' * jWindow)';
        end %for
    catch
        disperr_();
    end
    mrPc1 = (mrPc1') / nSamples;
    mrPc2 = (mrPc2') / nSamples;
end %func
