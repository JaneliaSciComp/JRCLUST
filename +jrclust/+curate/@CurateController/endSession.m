function res = endSession(obj)
    %ENDSESSION Finish curating and return results
    if nargout == 0
        success = obj.saveFiles();
        if ~success
            return;
        end
    end
    obj.isEnding = 1;
    obj.closeFigures();

    obj.currentSite = [];
    obj.selected = [];
    res = obj.cRes;
    obj.cRes = [];
    obj.isEnding = 0; % ended
end