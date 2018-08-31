%--------------------------------------------------------------------------
function mr = reshape_(vr, n1)
    % n1: leading dimension
    n = numel(vr);
    n2 = floor(n / n1);
    if n == (n1*n2)
        mr = reshape(vr, n1, n2);
    else
        mr = reshape(vr(1:(n1*n2)), n1, n2);
    end
end % function
