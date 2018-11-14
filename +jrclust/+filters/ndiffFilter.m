function samplesIn = ndiffFilter(samplesIn, nDiff_filt)
    %NDIFFFILTER
    if isempty(intersect(nDiff_filt, [1 2 3]))
        return;
    end

    if nDiff_filt == 1
        samplesIn(1:end-1,:) = diff(samplesIn);
        samplesIn(end, :) = 0;
    elseif nDiff_filt == 2
        samplesIn(2:end-2, :) = 2*diff(samplesIn(2:end-1, :)) + (samplesIn(4:end, :) - samplesIn(1:end-3, :));
        samplesIn([1, end-1, end], :) = 0;
    else % nDiff_filt == 3
        samplesIn(3:end-3, :) = 3*diff(samplesIn(3:end-2, :)) + 2*(samplesIn(5:end-1, :) - samplesIn(2:end-4, :)) + (samplesIn(6:end, :) - samplesIn(1:end-5, :));
        samplesIn([1, 2, end-2, end-1, end], :) = 0;
    end
end