function [success, errmsg] = closeFiltForWriting(obj)
    %CLOSEFILTFORWRITING Close the filtered file
    success = 0;
    errmsg = '';

    if obj.filteredFid < 3
        errmsg = 'Filtered file is not open';
        return;
    end

    success = (fclose(obj.filteredFid) == 0);
    obj.filteredFid = -1;
end

