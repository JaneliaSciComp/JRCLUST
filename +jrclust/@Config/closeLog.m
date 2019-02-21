function closeLog(obj)
    %CLOSELOG Close the log file
    if obj.logFid ~= -1
        try
            fclose(obj.logFid);
        catch ME
            warning('Failed to close log: %s', ME.message);
        end
    end

    obj.logFid = -1;
end

