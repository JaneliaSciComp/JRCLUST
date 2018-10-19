%--------------------------------------------------------------------------
function mn = uint2int_(mn)
    if isa(mn, 'uint16')
        mn = int16(single(mn)-2^15);
    elseif isa(mn, 'uint32')
        mn = int32(double(mn)-2^31);
    end
end
