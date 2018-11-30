function meanWf = zeroStart(meanWf)
    %ZEROSTART subtract the first row from all rows
    shape = size(meanWf);
    if numel(shape) ~= 2
        meanWf = reshape(meanWf, shape(1), []);
    end

    meanWf = bsxfun(@minus, meanWf, meanWf(1, :));
    if numel(shape) ~= 2
        meanWf = reshape(meanWf, shape);
    end
end