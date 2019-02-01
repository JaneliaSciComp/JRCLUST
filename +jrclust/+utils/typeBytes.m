function n = typeBytes(dataType)
    %TYPEBYTES Get number of bytes in dataType
    switch dataType
        case {'char', 'byte', 'int8', 'uint8'}
            n = 1;

        case {'int16', 'uint16'}
            n = 2;

        case {'single', 'float', 'int32', 'uint32'}
            n = 4;

        case {'double', 'int64', 'uint64'}
            n = 8;

        otherwise
            n = 0;
    end
end

