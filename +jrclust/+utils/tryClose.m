function tryClose(h)
    %TRYCLOSE Try to close a handle, fail gracefully
    try
        close(h);
    catch
    end
end
