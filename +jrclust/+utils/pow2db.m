function ydB = pow2db(y)
    ydB = (10.*log10(y) + 300) - 300;
end
