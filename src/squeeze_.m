function val = squeeze_(val, idimm)
    % val = squeeze_(val) : when squeezeing matrix, transpose if leading dimm is 1
    % val = squeeze_(val, idimm): permute specified dimension out
    s = size(val);

    if nargin >= 2
        dimm_ = [setdiff(1:ndims(val), idimm), idimm];
        val = permute(val, dimm_);
    elseif numel(s) == 2 && s(1) == 1
        val = val';
    else
        val = squeeze(val);
    end
end
