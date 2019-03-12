function [success, errmsg] = openFiltForWriting(obj)
    %OPENFILTFORWRITING Open filtered file for writing
    success = 0;
    errmsg = '';

    % first check we're not already open
    if obj.filteredFid > 2
        success = 1;
        return;
    end

    % next check we have the free space
    freeSpace = java.io.File(fileparts(obj.filtPath)).getFreeSpace();
    if freeSpace < obj.fSizeBytes % not enough space for another one
        errmsg = 'Not enough space on disk';
        return;
    end

    % still here? try to open the file for writing
    [obj.filteredFid, errmsg] = fopen(obj.filtPath, 'w');
    if obj.filteredFid == -1
        warning('Failed to open filtered file %s: %s', obj.filtPath, errmsg);
        return;
    end

    success = 1;
end

