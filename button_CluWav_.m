%--------------------------------------------------------------------------
function button_CluWav_(xyPos, vcButton)
    if strcmpi(vcButton, 'normal')
        event.Button = 1;
    elseif strcmpi(vcButton, 'alt')
        event.Button = 3;
    else
        return;
    end
    xPos = round(xyPos(1));
    S0 = get(0, 'UserData');

    switch(event.Button)
        case 1 %left click. copy clu and delete existing one
            S0 = update_cursor_(S0, xPos, 0);
        case 2 %middle, ignore
            return;
        case 3 %right click. paste clu
            S0 = update_cursor_(S0, xPos, 1);
    end

    figure_wait_(1);
    S0 = keyPressFcn_cell_(getCachedFig('FigWav'), {'t','c','i','v','e','f'}, S0); % 'z','j'

    set(0, 'UserData', S0);
    plot_raster_(S0);
    figure_wait_(0);
end %func
