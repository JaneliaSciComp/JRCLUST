function res = endSession(obj)
    %ENDSESSION Finish curating and return results
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    obj.isEnding = 1;
    if nargout == 0
        success = obj.saveFiles();
        if ~success
            obj.isEnding = 0; % cancelled
            return;
        end
    end
    obj.closeFigures();

    obj.currentSite = [];
    obj.selected = [];
    res = obj.cRes;
    obj.cRes = [];
    obj.isEnding = 0; % ended
end