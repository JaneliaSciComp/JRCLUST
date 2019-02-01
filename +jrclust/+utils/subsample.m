function [vals, vi] = subsample(vals, k, dim)
    %SUBSAMPLE Sample k items from vals, optionally along dim
    if nargin < 3
        dim = 2;
    end

    if isvector(vals)
        doReshape = isrow(vals);
        [vals, vi] = subsampleVec(vals(:), k);

        % reshape sampled values if necessary
        if doReshape
            vals = vals';
        end
    else
        if isempty(intersect(dim, 1:ndims(vals)))
            error('cannot sample along dimension %d', dim);
        end

        [vals, vi] = subsampleNd(vals, k, dim);
    end
end

%% LOCAL FUNCTIONS
function [vals, inds] = subsampleVec(vals, k)
    %SUBSAMPLEVEC Subsample `k` values from a vector
    if isempty(k) || k >= numel(vals)
        inds = 1:numel(vals);
    else
        inds = sort(randsample(1:numel(vals), k, 0));
    end

    vals = vals(inds);
end

function [vals, inds] = subsampleNd(vals, k, dim)
    %SUBSAMPLEND Subsample `k` rows (or cols, or ...) from an Nd array, N >= 2
    if isempty(k) || k >= size(vals, dim)
        inds = 1:size(vals, dim);
    else
        inds = sort(randsample(1:size(vals, dim), k, 0));
    end

    if dim == 1
        vals = vals(inds, :);
    else
        vals = vals(:, inds);
    end
end
