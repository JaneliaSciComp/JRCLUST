function [proj1, proj2, proj3] = pcProjectSpikes(sampledWindows, prVecs1, prVecs2, prVecs3)
    %PCPROJECTSPIKES Project selected spikes onto 1st, 2nd, 3rd principal components
    %   input:  sampledWindows is nSamples x nSpikes x nSites
    %   input:  prVecs{1,2,3} is nSamples x nSites
    %   output: proj{1,2,3} is nSites x nSpikes
    [nSamples, nSpikes, nSites] = size(sampledWindows);
    [proj1, proj2, proj3] = deal(zeros(nSpikes, nSites, 'single'));

    try
        for jSite = 1:nSites
            jWindow = jrclust.utils.meanSubtract(single(sampledWindows(:, :, jSite)));
            % project jWindow onto 1st principal vector for jSite
            proj1(:, jSite) = (prVecs1(:, jSite)' * jWindow)';
            if nargin > 2 && nargout > 1 % project jWindow onto 2nd principal vector for jSite
                proj2(:, jSite) = (prVecs2(:, jSite)' * jWindow)';
            end
            if nargin > 3 && nargout > 2 % project jWindow onto 3rd principal vector for jSite
                proj3(:, jSite) = (prVecs3(:, jSite)' * jWindow)';
            end
        end % for
    catch ME
        warning('error projecting spikes: %s', ME.message)
    end

    proj1 = (proj1') / nSamples; % why scale here?
    proj2 = (proj2') / nSamples;
    proj3 = (proj3') / nSamples;
end
