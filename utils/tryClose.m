%--------------------------------------------------------------------------
function tryClose(handle)
    % tryClose(handle) Try to close a handle and fail silently.
    try
        close(handle);
    catch
    end
end
