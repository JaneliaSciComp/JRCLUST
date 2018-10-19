%--------------------------------------------------------------------------
function hMsg = msgbox_open_(vcMessage)
    global fDebug_ui

    % if get_set_([], 'fDebug_ui', 0), hMsg = []; return; end
    if fDebug_ui==1, hMsg = []; return; end
    hFig = gcf;
    hMsg = msgbox(vcMessage);
    figure(hFig); %drawnow;
end
