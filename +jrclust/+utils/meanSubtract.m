function vals = meanSubtract(vals, dim, hFun)
    %MEANSUBTRACT Subtract mean (or other statistic) from vals
    if nargin < 2
        dim = 1;
    end
    if nargin < 3
        hFun = @mean;
    end

    if ~reallyisa(vals, 'single') && ~reallyisa(vals, 'double')
        vals = single(vals);
    end

    % have to reshape if ndims(vals) > 2
    shape = size(vals);
    if numel(shape) > 2
        vals = reshape(vals, shape(1), []);
    end

    % compute mean and subtract it from vals
    vals = bsxfun(@minus, vals, hFun(vals, dim));

    % restore original shape
    if numel(shape) > 2
        vals = reshape(vals, shape);
    end
end

%% LOCAL FUNCTIONS
function flag = reallyisa(val, clsname)
    %REALLYISA Like isa, but checks the underlying class if a gpuArray
    try
        if isa(val, 'gpuArray')
            flag = strcmpi(clsname, classUnderlying(val));
        else
            flag = isa(val, clsname);
        end
    catch
        flag = 0;
    end
end
