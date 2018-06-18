%--------------------------------------------------------------------------
function mrPc1 = pc1_tr_(tn)
    % returns first principal component across sites

    mr0 = single(squeeze_(tn(:,1,:)));
    vrPv0 = zscore_(pca(mr0', 'NumComponents', 1));
    dimm_tn = size(tn);
    mr = single(reshape(tn, dimm_tn(1), []));
    mr = bsxfun(@minus, mr, mean(mr));
    mrPc1 = reshape(vrPv0' * mr, dimm_tn(2:3));
end %func
