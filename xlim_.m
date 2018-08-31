%--------------------------------------------------------------------------
function xlim_(arg1, arg2)
    % ylim function
    % ylim_(lim_)
    % ylim_(hAx, lim_)
    if nargin==1
        [hAx_, lim_] = deal(gca, arg1);
    else
        [hAx_, lim_] = deal(arg1, arg2);
    end
    if any(isnan(lim_)), return; end
    try
        xlim(hAx_, sort(lim_));
    catch
        disperr_();
    end
end % function
