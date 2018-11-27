function [vals, vi] = subsample(vals, k, dim)
    %SUBSAMPLE Sample k items from vals, optionally along dim
    shape = size(vals);

    if ~any(shape == 1) || nargin == 3
        [vals, vi] = subsampleMat(vals, k, dim);
    else
        [vals, vi] = subsampleVec(vals, k);
    end
end

%% LOCAL FUNCTIONS
function [vals, vi] = subsampleVec(vals, k)
    %SUBSAMPLEVEC
    vi = [];

    if numel(vals) > k
        % nSkip = floor(numel(vals)/n);
        % vi = 1:nSkip:numel(vals);
        vi = sort(randsample(1:numel(vals), k, false));
        vals = vals(vi);
    end
end

function [vals, vi] = subsampleMat(vals, n, dim)
    %SUBSAMPLEMAT
    if nargin < 3
        dim = 2;
    end

    if isempty(intersect(dim, [1 2]))
        error('bad dimension');
    end

    if isempty(n)
        return;
    end

    n = size(vals, dim);
    nSkip = max(floor(n / n), 1); %?
    vi = 1:nSkip:n;

    if nSkip == 1
        return;
    end

    vi = vi(1:n);

    if dim == 1
        vals = vals(vi, :);
    else
        vals = vals(:, vi);
    end

    if nargout>=2
        if n > n
            vi = 1:nSkip:n;
            vi = vi(1:n);
        else
            vi = 1:n;
        end
    end
end
