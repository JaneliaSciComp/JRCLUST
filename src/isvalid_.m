%--------------------------------------------------------------------------
function flag = isvalid_(h)
    if isempty(h), flag = 0; return ;end
    try
        flag = isvalid(h);
    catch
        flag = 0;
    end
end %func
