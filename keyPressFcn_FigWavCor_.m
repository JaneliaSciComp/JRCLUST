%--------------------------------------------------------------------------
function S0 = keyPressFcn_FigWavCor_(~, event)
    S0 = get(0, 'UserData');
    switch lower(event.Key)
        case 'm' %merge
            S0 = ui_merge_(S0);
        case 's' %split
            S0 = auto_split_(1, S0); %multi
        case {'d', 'backspace', 'delete'} %delete
            S0 = ui_delete_(S0);
    end %switch

    if nargout == 0
        set(0, 'UserData', S0);
    end
end %func
