function openFilt(obj)
    %OPENFILT Open the filtered recording file for reading
    if obj.filtIsOpen
        return;
    end

    dshape = [obj.hCfg.nSites, obj.nSamples];
    obj.filteredData = memmapfile(obj.filtPath, 'Offset', 0, ...
                                  'Format', {obj.dataType, dshape, 'Data'}, ...
                                  'Writable', false);
    obj.filtIsOpen = 1;
end