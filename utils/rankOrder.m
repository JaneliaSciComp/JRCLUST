%--------------------------------------------------------------------------
function [ranking, indices] = rankOrder(vals, ordering)
    % Return index of each element of VALS, if VALS were a sorted array; an unsorter.
    % i.e., if sortedVals == vals(indices), then vals == sortedVals(ranking)
    % warning: 32 bit addressing (JJJ)

    if nargin < 2
        ordering = 'ascend';
    end

    n = numel(vals);

    [~, indices] = sort(vals, ordering);

    if isGpu_(vals)
        ranking = zeros(n, 1, 'int32', 'gpuArray');
        ranking(indices) = 1:n;
    else
        ranking = zeros(n, 1, 'int32');
        ranking(indices) = 1:n;
    end
end % function
