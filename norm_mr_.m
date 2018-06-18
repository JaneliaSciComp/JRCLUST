%--------------------------------------------------------------------------
function [mr, vr] = norm_mr_(mr)
    % Normalize each columns
    vr = sqrt(sum(mr .^ 2));
    mr = bsxfun(@rdivide, mr, vr);
end %func
