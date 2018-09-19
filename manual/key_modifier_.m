%--------------------------------------------------------------------------
function flag = key_modifier_(event, vcKey)
    % Check for shift, alt, ctrl press
    try
        flag = any(strcmpi(event.Modifier, vcKey));
    catch
        flag = 0;
    end
end %func
