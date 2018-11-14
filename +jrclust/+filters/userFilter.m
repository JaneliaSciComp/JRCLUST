function samplesIn = userFilter(samplesIn, filtKernel)
    %USERFILTER Convolve samplesIn with user's custom kernel

    if ~isempty(filtKernel)
        for i = 1:size(samplesIn, 2)
            samplesIn(:, i) = conv(samplesIn(:, i), filtKernel, 'same');
        end
    end
end

