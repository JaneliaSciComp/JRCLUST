%--------------------------------------------------------------------------
function gui_test_(P, vcFig, csMenu_skip)
    if nargin<3, csMenu_skip = {}; end
    drawnow;
    hFig = get_fig_(vcFig);
    keyPress_fig_(hFig, get_keyPress_('all'));

    % Menu test
    menu_test_(hFig, csMenu_skip);

    try
        close(hFig); %close traces figure. other figures may remain
    catch
        ;
    end
end %func
