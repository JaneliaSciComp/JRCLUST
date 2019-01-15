function flag = keyMod(hEvent, kMod)
    %KEYMOD Check for a key modifier
    try
        flag = any(strcmpi(hEvent.Modifier, kMod));
    catch
        flag = 0;
    end
end
