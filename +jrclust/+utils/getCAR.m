function car = getCAR(samplesIn, CARMode, ignoreSites)
    %GETCAR Get common average reference for samples (mean or median)
    if ~isempty(ignoreSites)
        samplesIn(:, ignoreSites) = [];
    end

    if strcmpi(CARMode, 'median')
        car = median(samplesIn, 2);
    else
        car = mean(samplesIn, 2);
    end
end
