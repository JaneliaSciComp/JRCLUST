function val = inProgress(obj)
    %INPROGRESS Return true if in progress (not finished or errored)
    val = ~(obj.isCompleted || obj.isError) && ~isempty(obj.hCfg);
end