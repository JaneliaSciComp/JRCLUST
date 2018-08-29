%--------------------------------------------------------------------------
function vr = mr2vr_sub2ind_(mr, vi1, vi2)
    if isempty(mr), vr = []; return; end
    if isempty(vi1), vi1 = 1:size(mr,1); end
    if isempty(vi2), vi2 = 1:size(mr,2); end
    vr = mr(sub2ind(size(mr), vi1(:), vi2(:)));
end % function
