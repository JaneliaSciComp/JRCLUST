function flag = isvalid(h)
    %ISVALID Like isvalid, but checks for empty first
    if isempty(h)
        flag = 0;
        return;
    end

    try
        flag = isvalid(h);
    catch
        flag = 0;
    end
end
