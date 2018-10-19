%--------------------------------------------------------------------------
function [icl, x, y] = log_ztran_(x, y, x_cut, y_cut)
    % [icl, x, y] = detrend_ztran_(x, y, x_cut, y_cut)
    % [icl, x, y] = detrend_ztran_(x, y, n_icl)
    if nargin == 3
        n_icl = x_cut;
    else
        n_icl = [];
    end

    x = log10(x(:));
    y = log10(y(:));

    vlValid = isfinite(x) & isfinite(y);
    y(vlValid) = zscore_(y(vlValid));

    if isempty(n_icl)
        icl = find(x>=x_cut & y>=y_cut);
    else
        [~, ix] = sort(y(vlValid), 'descend');
        viValid = find(vlValid);
        icl = viValid(ix(1:n_icl));
    end

    if nargout==0
        figure; plot(x,y,'.', x(icl),y(icl),'ro'); grid on;
    end
end %func
