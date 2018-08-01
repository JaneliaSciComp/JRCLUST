%--------------------------------------------------------------------------
function keyPressFcn_FigWavCor_(hObject, event)
    S0 = get(0, 'UserData');
    switch lower(event.Key)
        case 'm' %merge
        ui_merge_(S0);
        case 's' %split
        auto_split_(1); %multi
        case {'d', 'backspace', 'delete'} %delete
        ui_delete_(S0);
    end %switch
end %func
