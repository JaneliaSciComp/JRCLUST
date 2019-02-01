% function [vlKeep_ref, channelMeansMAD] = carReject(vrWav_mean1, P)
function [keepMe, channelMeansMAD] = carReject(channelMeans, blankPerMs, blankThresh, sampleRate)
    %CARREJECT

    channelMeansMAD = [];

    channelMeans = single(channelMeans);
    blankWindow = ceil(sampleRate*blankPerMs/1000);
    
    if blankWindow > 0
        if nargout < 2
            keepMe = threshMAD(abs(channelMeans), blankThresh);
        else
            [keepMe, channelMeansMAD] = threshMAD(abs(channelMeans), blankThresh);
        end

        if blankWindow > 1 % clear out neighbors of blanked samples
            overThresh = find(~keepMe);
            for crossing = overThresh'
                left = max(crossing - ceil(blankWindow/2), 1);
                right = min(crossing + ceil(blankWindow/2), numel(keepMe));
                keepMe(left:right) = 0;
            end
        end
%     else
%         channelMeansBinned = std(padReshape(channelMeans, blankWindow), 1, 1);
%         if nargout < 2
%             keepMe = threshMAD(channelMeansBinned, blankThresh);
%         else
%             [keepMe, channelMeansMAD] = threshMAD(channelMeansBinned, blankThresh);
%             channelMeansMAD = expand_vr_(channelMeansMAD, blankWindow, size(channelMeans));
%         end
% 
%         keepMe = expand_vr_(keepMe, blankWindow, size(channelMeans));
    else
        [keepMe, channelMeansMAD] = deal([]);
    end
end

%% LOCAL FUNCTIONS
function [keepMe, channelMeans] = threshMAD(channelMeans, madThresh)
    %THRESHMAD Find which channel means fall within MAD threshold
    nSubsamples = 300000;

    % estimate MAD of channelMeans
    med = median(jrclust.utils.subsample(channelMeans, nSubsamples));
    channelMeans = channelMeans - med; % deviation from the median
    mad = median(abs(jrclust.utils.subsample(channelMeans, nSubsamples)));

    if isempty(madThresh) || madThresh == 0
        keepMe = true(size(channelMeans));
    else
        keepMe = channelMeans < mad*madThresh;
    end

    if nargout >= 2
        channelMeans = channelMeans/mad; % MAD unit
    end
end

% function mat = padReshape(vec, nwin)
%     %PADRESHAPE Reshape a vector to a matrix, padding if necessary
%     nbins = ceil(numel(vec)/nwin);
%     vec(nbins*nwin) = 0; % pad the end with zeros
%     mat = reshape(vec(1:nbins*nwin), nwin, nbins);
% end
% 
% function vr1 = expand_vr_(vr, nwin, shape)
%     if islogical(vr)
%         vr1 = false(shape);
%     else
%         vr1 = zeros(shape, 'like', vr);
%     end
%     vr = repmat(vr(:)', [nwin, 1]);
%     vr = vr(:);
%     [n,n1] = deal(numel(vr), numel(vr1));
%     if n1 > n
%         vr1(1:n) = vr;
%         vr1(n+1:end) = vr1(n);
%     elseif numel(vr1) < n
%         vr1 = vr(1:n1);
%     else
%         vr1 = vr;
%     end
% end

