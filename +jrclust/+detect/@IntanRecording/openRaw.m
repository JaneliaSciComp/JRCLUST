function openRaw(obj)
    %OPENRAW Open the raw recording file for reading
    if obj.rawIsOpen
        return;
    end

    obj.rawFid = fopen(obj.rawPath, 'r');
    fseek(obj.rawFid, obj.headerOffset, 'bof');
    obj.rawIsOpen = 1;
end