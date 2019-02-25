function openRaw(obj)
    %OPENRAW Open the raw recording file for reading
    if obj.rawIsOpen
        return;
    end

    obj.rawData = memmapfile(obj.rawPath, 'Offset', obj.headerOffset, ...
                             'Format', {obj.dataType, obj.dshape, 'Data'}, ...
                             'Writable', false);
    obj.rawIsOpen = 1;
end