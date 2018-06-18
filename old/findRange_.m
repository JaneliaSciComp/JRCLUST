%--------------------------------------------------------------------------
% 10/11/17 JJJ: new and faster. tested
function [vii, ic] = findRange_(vi, a, b, i, n)
    % a: start value
    % b: end value
    % i: index of vi to start searching
    % vii: index of vi that is between a and b
    % vii = [];
    ic = 1;
    ia = 0;
    vii = [];
    if b < vi(1), return; end % exception condition added
    if a > vi(end), ic = n; return; end
    while 1
        v_ = vi(i);
        if v_ < a, ic = i; end
        if ia==0
            if v_ >= a, ia = i; end
        else
            if v_ > b, ib = i-1; break; end
        end
        if i<n
            i = i + 1;
        else
            ib = n;
            break;
        end
    end
    if ia > 0
        vii = ia:ib;
    end
end %func
