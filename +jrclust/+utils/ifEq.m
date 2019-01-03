function out = ifEq(pred, onTrue, onFalse)
    %IFEQ If pred is true, return onTrue, otherwise return onFalse
    if (pred)
        out = onTrue;
    else
        out = onFalse;
    end
end
