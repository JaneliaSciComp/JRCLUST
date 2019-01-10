function hMsgbox = qMsgBox(msg, fBlock, fModal)
    %QMSGBOX msgbox with some reasonable defaults
    hMsgbox = [];

    if nargin < 2
        fBlock = false;
    end
    if nargin < 3
        fModal = false;
    end

    if fBlock
        uiwait(msgbox(msg, 'modal'));
    else
        try
            hMsgbox = msgbox(msg, jrclust.utils.ifEq(fModal, 'modal', 'non-modal'));
        catch
        end
    end
end
