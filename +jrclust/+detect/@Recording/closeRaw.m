function closeRaw(obj)
    %CLOSERAW Close the raw recording file, clear its data
    if ~obj.rawIsOpen
        return;
    end
    obj.rawData = [];
    obj.rawIsOpen = 0;
end