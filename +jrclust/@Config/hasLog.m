function flag = hasLog(obj, logTag)
    %HASLOG Do we have a log entry with this tag?
    flag = ischar(logTag) && isKey(obj.logEntries, logTag);
end