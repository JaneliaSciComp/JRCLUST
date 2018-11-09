function ok = ismatrixnum(val)
    %ISMATRIXNUM Return true iff a val is a numeric matrix
    ok = isnumeric(val) && ismatrix(val);
end

