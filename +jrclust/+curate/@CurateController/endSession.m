function res = endSession(obj)
    %ENDSESSION Finish curating and return results
    if obj.isWorking
        dlgAns = questdlg('An operation is in progress. Really quit?', 'Operation in progress', 'No');
        if ~strcmp(dlgAns, 'Yes')
            return;
        end
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