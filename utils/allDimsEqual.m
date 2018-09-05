%--------------------------------------------------------------------------
function isequal = allDimsEqual(val, dims)
    valSize = size(val);
    nDims = numel(dims);

    isequal = (numel(valSize) == nDims) && all(valSize == dims);
end
