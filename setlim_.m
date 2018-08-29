%--------------------------------------------------------------------------
% 17/9/13 JJJ: Created and tested
function vr = setlim_(vr, lim_)
    % Set low and high limits
    vr = min(max(vr, lim_(1)), lim_(2));
end % function
