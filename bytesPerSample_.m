%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function n = bytesPerSample_(vcDataType)
    % Return number of bytes per data type

    switch lower(vcDataType)
        case {'char', 'byte', 'int8', 'uint8'}
        n = 1;
        case {'int16', 'uint16'}
        n = 2;
        case {'single', 'float', 'int32', 'uint32'}
        n = 4;
        case {'double', 'int64', 'uint64'}
        n = 8;
        otherwise
        n = [];
        fprintf(2, 'Unsupported data type: %s\n', vcDataType);
    end
end %func
