%--------------------------------------------------------------------------
function flag = equal_vr_(vr1, vr2)
    % return true if and only if sizes are equal and all elements are equal
    flag = (all(size(vr1) == size(vr2)) && all(vr1(:) == vr2(:)));
    % if all(size(vr1) == size(vr2))
    %     ml = vr1 == vr2;
    %     flag = all(ml(:));
    % else
    %     flag = 0;
    % end
end %func
