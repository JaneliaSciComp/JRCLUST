%--------------------------------------------------------------------------
function tr = meanSubt_tr_(tr)
    dimm = size(tr);
    mr = reshape(tr,size(tr,1),[]);
    mr = bsxfun(@minus, mr, mean(mr));
    tr = reshape(mr, dimm);
end %func
