%--------------------------------------------------------------------------
function flag = keyModifier(event, vcKey)
    % Check for shift, alt, ctrl press
    try
        flag = any(strcmpi(event.Modifier, vcKey));
    catch
        flag = 0;
    end
end %func
