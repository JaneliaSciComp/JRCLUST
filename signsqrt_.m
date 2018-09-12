%--------------------------------------------------------------------------
% Return sign-preserving square root
% 7/26/17 JJJ: Code cleanup and test
function x = signsqrt_(x)
    x = sign(x) .* sqrt(abs(x));
end %func
