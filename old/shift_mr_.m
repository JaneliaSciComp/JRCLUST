%--------------------------------------------------------------------------
function mr = shift_mr_(mr, n)
    if n<0
        n=-n;
        mr(1:end-n, :) = mr(n+1:end, :);
    elseif n>0
        mr(n+1:end, :) = mr(1:end-n, :);
    else %n==0
        return;
    end
end %func
