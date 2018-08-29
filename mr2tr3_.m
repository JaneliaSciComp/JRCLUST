%--------------------------------------------------------------------------
function [tr, miRange] = mr2tr3_(mr, spkLim, viTime, viSite, fMeanSubt)
    % tr: nSamples x nSpikes x nChans

    if nargin < 4
        viSite = [];
    end %faster indexing
    if nargin < 5
        fMeanSubt = 0;
    end

    % JJJ 2015 Dec 24
    % vr2mr2: quick version and doesn't kill index out of range
    % assumes vi is within range and tolerates spkLim part of being outside
    % works for any datatype
    if isempty(viTime), tr=[]; return; end
    [N, M] = size(mr);
    if ~isempty(viSite), M = numel(viSite); end
    if iscolumn(viTime), viTime = viTime'; end

    viTime0 = [spkLim(1):spkLim(end)]'; %column
    miRange = bsxfun(@plus, int32(viTime0), int32(viTime));
    miRange = min(max(miRange, 1), N);
    miRange = miRange(:);

    if isempty(viSite)
        tr = mr(miRange,:);
    else
        tr = mr(miRange, viSite);
    end
    tr = reshape(tr, [numel(viTime0), numel(viTime), M]);

    if fMeanSubt
        %     trWav1 = single(permute(trWav1, [1,3,2]));
        tr = single(tr);
        dimm1 = size(tr);
        tr = reshape(tr, size(tr,1), []);
        tr = bsxfun(@minus, tr, mean(tr)); %mean subtract
        tr = reshape(tr, dimm1);
    end
end % function
