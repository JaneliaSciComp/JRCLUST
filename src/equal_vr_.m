%--------------------------------------------------------------------------
function flag = equal_vr_(vr1, vr2)
    if all(size(vr1) == size(vr2))
        ml = vr1 == vr2;
        flag = all(ml(:));
    else
        flag = 0;
    end
end %func
