%--------------------------------------------------------------------------
function [icl, x, z] = detrend_ztran_(x0, y0, x_cut, y_cut)
    % [icl, x, y] = detrend_ztran_(x, y, x_cut, y_cut)
    % [icl, x, y] = detrend_ztran_(x, y, n_icl)
    if nargin == 3
        n_icl = x_cut;
    else
        n_icl = [];
    end
    y_thresh_mad = 20;
    x = log10(x0(:));
    y = (y0(:));
    vlValid_x = x>=x_cut & x <= -1;
    z_max = 100;

    % compute y_mad
    y_med = median(y(vlValid_x));
    y_mad = y - y_med;
    y_mad = y_mad / median(y_mad);
    vlValid_y = abs(y_mad) < y_thresh_mad;
    % vl_y0 = y0<=0;
    % y = y0;
    % y(vl_y0) = min(y(~vl_y0));
    % y = log10(y0(:));
    % max_y = max(y);

    vlValid = isfinite(x) & isfinite(y);
    % viDetrend = find(vlValid  & (y < max_y/4) & (x >= x_cut));
    viDetrend = find(vlValid  & vlValid_y & vlValid_x);
    x1 = x(viDetrend);
    y1 = y(viDetrend);
    xx1 = [x1+eps('single'), ones(numel(x1),1)];
    m = xx1 \ y1; %determine slope based on acceptable units

    xx = [x+eps('single'), ones(numel(x),1)];
    y2 = y - xx * m; %detrend data
    mu2 = mean(y2(viDetrend));
    sd2 = std(y2(viDetrend));
    z = (y2 - mu2) / sd2; %z transformation
    z(z>z_max) = z_max;

    if isempty(n_icl)
        icl = find(x>=x_cut & z>=10^y_cut);
    else
        [~, ix] = sort(z(vlValid), 'descend');
        viValid = find(vlValid);
        icl = viValid(ix(1:n_icl));
    end

    if nargout==0
        figure; plot(x,z,'.', x(icl),z(icl),'ro'); grid on;
        title(sprintf('%d clu', numel(icl)));
    end
end %func
