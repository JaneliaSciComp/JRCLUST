function [multiBounds, multiRange, multiEdges] = sample_skip_(windowBounds, nSamplesTotal, nTimeTraces)
    if nTimeTraces == 1 || isempty(nTimeTraces)
        multiBounds = {windowBounds};
        multiRange = windowBounds(1):windowBounds(end);
        multiEdges = [];
        return;
    end

    nSkip = floor(nSamplesTotal / nTimeTraces);
    multiBounds = arrayfun(@(i) windowBounds + (i-1)*nSkip, 1:nTimeTraces, 'UniformOutput', 0);

    for i = 1:nTimeTraces
        lim1 = mod(multiBounds{i} - 1, nSamplesTotal) + 1;
        if lim1(1) > lim1(2) % lim1 runs out of bounds
            lim1 = [1, diff(windowBounds) + 1]; % is this correct?
        end
        multiBounds{i} = lim1;
    end

    if nargout >= 2
        multiRange = jrclust.utils.neCell2mat(cellfun(@(x) x(1):x(2), multiBounds, 'UniformOutput', 0));
    end

    if nargout >= 3 %compute the number of samples
        multiEdges = cumsum(cellfun(@(x) diff(x) + 1, multiBounds));
        multiEdges = [1, multiEdges(1:end-1)];
    end
end
