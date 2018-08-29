%--------------------------------------------------------------------------
function hMsgbox = msgbox_(csMsg, fBlock, fModal)
    % msgbox. Don't display if fDebug_ui is set
    hMsgbox = [];
    if nargin<2, fBlock = 0; end
    if nargin<3, fModal = 0; end
    global fDebug_ui
    if fDebug_ui==1, return; end
    if fBlock
        uiwait(msgbox(csMsg, 'modal'));
    else
        try
            hMsgbox = msgbox(csMsg, ifeq_(fModal, 'modal', 'non-modal'));
        catch
            ;
        end
    end
end % function
