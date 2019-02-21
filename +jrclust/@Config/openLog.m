function openLog(obj, logFilename)
    %OPENLOG Open the log file
    if nargin < 2 || isempty(logFilename)
        logFilename = [datestr(now(), 30) '.log'];
    end

    % no directory specified in filename
    if isempty(fileparts(logFilename))
        logDirname = fullfile(obj.outputDir, '.jrc');
        if exist(logDirname, 'dir') ~= 7
            try
                mkdir(logDirname);
            catch ME
                warning('Could not create log directory %s: %s', logDirname, ME.message);
                obj.logFid = -1;
                return;
            end
        end

        logFilename = fullfile(logDirname, logFilename);
    end
    
    [obj.logFid, errmsg] = fopen(logFilename, 'w');
    if obj.logFid == -1
        warning('Failed to open log: %s', errmsg);
    end
end

