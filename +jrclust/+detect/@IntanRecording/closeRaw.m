function closeRaw(obj)
    %CLOSERAW Close the raw recording file, clear its data
    if obj.rawFid == -1
        return;
    end

    try
        fclose(obj.rawFid);
    catch ME
    end
    obj.rawIsOpen = 0;
end