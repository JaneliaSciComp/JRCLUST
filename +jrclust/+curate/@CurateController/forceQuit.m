function cRes = forceQuit(obj)
    %FORCEQUIT Like endSession, but ignores isWorking
    if obj.isWorking
        obj.isWorking = 0;
    end

    if nargout > 0
        cRes = obj.endSession();
    else
        obj.endSession();
    end
end

