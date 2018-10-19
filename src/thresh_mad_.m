%--------------------------------------------------------------------------
function [vl, vr] = thresh_mad_(vr, thresh_mad)
    % single sided, no absolute value

    nsubs = 300000;
    offset = median(subsample_vr_(vr, nsubs));
    vr = vr - offset; %center the mean
    factor = median(abs(subsample_vr_(vr, nsubs)));
    if isempty(thresh_mad) || thresh_mad==0
        vl = true(size(vr));
    else
        vl = vr < factor * thresh_mad;
    end
    if nargout>=2
        vr = vr / factor; %MAD unit
    end
end %func
