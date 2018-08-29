%--------------------------------------------------------------------------
function [y_, z] = detrend_(x, y, vi, fQuadratic)
    % perform quadratic detrend and return z score
    if nargin<4, fQuadratic = 1; end
    if nargin<3, vi = 1:numel(x); end
    x=x(:);
    y=y(:);
    % y = log10(y);
    eps_ = eps(class(x));
    if fQuadratic
        xx = [x.^2+eps_, x + eps_, ones(numel(x),1)];
    else %linear detrending
        xx = [x + eps_, ones(numel(x),1)];
    end
    m = xx(vi,:) \ y(vi);
    y_ = y - xx * m;
    if nargout>=2
        z = (y_ - nanmean(y_(vi))) / nanstd(y_(vi));
    end
end % function
