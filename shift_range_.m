%--------------------------------------------------------------------------
function [cvi1, cvi2] = shift_range_(nT, nShift, viShift)
    % return ranges of two matrix to be time shifted
    % [cvi1, cvi2] = shift_range_(nT, nShift): -nShift:nShift
    % [cvi1, cvi2] = shift_range_(nT, vnShift)
    if ~isempty(nShift)
        [cvi1, cvi2] = deal(cell(nShift*2+1, 1));
        viShift = -nShift:nShift;
    else
        [cvi1, cvi2] = deal(cell(numel(viShift), 1));
    end
    viRange = 1:nT;
    for iShift_ = 1:numel(viShift)
        iShift = viShift(iShift_);
        iShift1 = -round(iShift/2);
        iShift2 = iShift + iShift1;
        viRange1 = viRange + iShift1;
        viRange2 = viRange + iShift2;
        vl12 = (viRange1>=1 & viRange1<=nT) & (viRange2>=1 & viRange2<=nT);
        cvi1{iShift_} = viRange1(vl12);
        cvi2{iShift_} = viRange2(vl12);
    end
end % function
