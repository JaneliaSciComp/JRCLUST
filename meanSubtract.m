%--------------------------------------------------------------------------
% 8/6/17 JJJ: Tested and documented
function traces = meanSubtract(traces, dim, hFunc)
    % subtract mean for mr or tr
    if nargin < 2
        dim = 1;
    end
    if nargin < 3
        hFunc = @mean;
    end

    if ~isa_(traces, 'single') && ~isa_(traces, 'double')
        traces = single(traces);
    end

    traceDims = size(traces);
    if numel(traceDims) > 2
        traces = reshape(traces, traceDims(1), []);
    end

    traces = bsxfun(@minus, traces, hFunc(traces, dim));
    if numel(traceDims) > 2
        traces = reshape(traces, traceDims);
    end
end % function
