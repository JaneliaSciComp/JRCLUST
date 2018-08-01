%--------------------------------------------------------------------------
function flag = tryIsValid(h)
    if isempty(h)
        flag = 0;
    else
        try
            flag = isvalid(h);
        catch
            flag = 0;
        end
    end
end % func
