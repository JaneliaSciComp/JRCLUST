function z = nanlog10(z)
    %NANLOG10 Replace negative values with nan to force real log10
    z(z <= 0) = nan;
    z = log10(z);
end