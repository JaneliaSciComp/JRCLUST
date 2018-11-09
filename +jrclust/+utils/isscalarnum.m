function ok = isscalarnum(val)
    %ISSCALARNUM Return true iff a val is a numeric scalar
    ok = isnumeric(val) && isscalar(val);
end

