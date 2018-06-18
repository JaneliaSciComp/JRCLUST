%--------------------------------------------------------------------------
function limit1 = limit_(vn)
    vn = vn(:);
    limit1 = [min(vn), max(vn)];
end
