function success = writeFilt(obj, filteredData)
    %WRITEFILT Write filtered data to a file
    success = 0;

    if ~obj.openFiltForWriting()
        warning('Failed to open filtered file');
        return;
    end

    count = fwrite(obj.filteredFid, filteredData, ['*' obj.dataType]);
    success = (count == numel(filteredData));
end

