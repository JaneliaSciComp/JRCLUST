%--------------------------------------------------------------------------
% Return sign-preserving log
% 7/26/17 JJJ: Code cleanup and test
function x = signlog_(x, e)
    if nargin<2, e = .0001; end
    x = sign(x) .* log(abs(x + e));
end %func
