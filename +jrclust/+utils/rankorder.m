function [ranking, argsort] = rankorder(vals, direction)
    %RANKORDER Rank vals in order of sorting direction.
    % Returns ranking and permutation needed to order vals
    if nargin < 2
        direction = 'ascend';
    end

    n = numel(vals);
    [~, argsort] = sort(vals, direction);

    ranking = zeros(n, 1, 'int32');
    ranking(argsort) = 1:n;
end
