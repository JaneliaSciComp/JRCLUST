%--------------------------------------------------------------------------
function mr = vr2mr2_(vr, viRow, spkLim, viCol)
    if nargin<4, viCol = []; end
    % JJJ 2015 Dec 24
    % vr2mr2: quick version and doesn't kill index out of range
    % assumes vi is within range and tolerates spkLim part of being outside
    % works for any datatype

    % prepare indices
    if size(viRow,2)==1, viRow=viRow'; end %row
    viSpk = int32(spkLim(1):spkLim(end))';
    miRange = bsxfun(@plus, viSpk, int32(viRow));
    miRange = min(max(miRange, 1), numel(vr));
    if isempty(viCol)
        mr = vr(miRange); %2x faster
    else
        mr = vr(miRange, viCol); %matrix passed to save time
    end
end %func
