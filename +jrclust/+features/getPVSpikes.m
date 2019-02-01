function [prVecs1, prVecs2, prVecs3] = getPVSpikes(sampledWindows)
    %GETPVSPIKES Get first 3 principal vectors for sampledWindows on each site
    %   input:  sampledWindows is nSamples x nSpikes x nSites
    %   output: prVecs{1,2,3} is nSamples x nSites
%     mrPv_global = get_(S0, 'mrPv_global');
%     if ~isempty(mrPv_global)
%         % show global pca
%         mrPv1 = repmat(mrPv_global(:,1), [1, nSites]);
%         mrPv2 = repmat(mrPv_global(:,2), [1, nSites]);
%     else

    if ismatrix(sampledWindows)
        nSamples = size(sampledWindows, 1);
        nSites = 1;
    else
        [nSamples, ~, nSites] = size(sampledWindows);
    end

    [prVecs1, prVecs2, prVecs3] = deal(zeros(nSamples, nSites, 'single'));
    for iSite = 1:nSites
        prVecs = jrclust.features.getPVSamples(sampledWindows(:, :, iSite));
        prVecs1(:, iSite) = prVecs(:, 1);
        prVecs2(:, iSite) = prVecs(:, 2);
        prVecs3(:, iSite) = prVecs(:, 3);
    end
%     end
end
