%--------------------------------------------------------------------------
function flag = isVisible_(hObj)
    flag = strcmpi(get(hObj, 'Visible'), 'on');
end %func
