function vals = u2i(vals)
    %U2I Convert an unsigned integer to ... something
    if isa(vals, 'uint16')
        vals = int16(single(vals)-2^15);
    elseif isa(vals, 'uint32')
        vals = int32(double(vals)-2^31);
    end
end
