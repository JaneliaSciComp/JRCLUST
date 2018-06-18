%--------------------------------------------------------------------------
function axis_(arg1, arg2)
    if nargin==1
        [hAx_, lim_] = deal(gca, arg1);
    else
        [hAx_, lim_] = deal(arg1, arg2);
    end
    if ischar(lim_)
        axis(hAx_, lim_);
        return;
    end
    if any(isnan(lim_)), return; end
    try
        axis(hAx_, [sort(lim_(1:2)), sort(lim_(3:4))]);
    catch
        disperr_();
    end
end %func
