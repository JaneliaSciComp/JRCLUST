function updateLog(obj, logTag, message, start, finish)
    %UPDATELOG Save and optionally output a message
    if nargin < 4
        start = 0;
    end
    if nargin < 5
        finish = 0;
    end

    if obj.hasLog(logTag)
        logEntry = obj.logEntries(logTag);
        logEntry.lastUpdate = now();
    else
        openedOn = now();
        logEntry = struct('message', message, ...
                          'opened', openedOn, ...
                          'lastUpdate', openedOn);
    end

    printMsg = sprintf('%s %s', datestr(now(), 31), message);
    if start
        logEntry.tic = tic();
        printMsg = sprintf('%s...', printMsg);
    end
    if finish && isfield(logEntry, 'tic')
        logEntry.toc = toc(logEntry.tic);
        printMsg = sprintf('%s (took %0.2f s)', printMsg, logEntry.toc);
    end

    obj.logEntries(logTag) = logEntry;

    if obj.verbose
        fprintf('%s\n', printMsg);
    end

    if obj.logFid > -1
        fprintf(obj.logFid, '%s\n', printMsg);
    end
end