function rerun(obj)
    %RERUN Rerun commands
    if obj.isError
        error(obj.errMsg);
    else
        obj.isCompleted = 0;
        obj.run();
    end
end