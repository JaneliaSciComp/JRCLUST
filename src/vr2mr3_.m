%--------------------------------------------------------------------------
function [mr, miRange] = vr2mr3_(vr, vi, spkLim)
    % JJJ 2015 Dec 24
    % vr2mr2: quick version and doesn't kill index out of range
    % assumes vi is within range and tolerates spkLim part of being outside
    % works for any datatype

    % prepare indices
    if size(vi,2)==1, vi=vi'; end %row
    viSpk = int32(spkLim(1):spkLim(end))';
    miRange = bsxfun(@plus, viSpk, int32(vi));
    miRange(miRange<1) = 1;
    miRange(miRange > numel(vr)) = numel(vr); %keep # sites consistent
    % miRange = int32(miRange);

    % build spike table
    nSpks = numel(vi);
    mr = reshape(vr(miRange(:)), [], nSpks);
end %func
