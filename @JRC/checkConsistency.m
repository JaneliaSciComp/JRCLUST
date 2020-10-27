function success = checkConsistency(obj)
%CHECKCONSISTENCY Check consistency of hClust and recover if necessary.
%   Warn user if data is found to be in an incconsistent state. Try to
%   recover the data automatically.
success = 1;

if ~isfield(obj.res, 'hClust')
    return;
end

if ~isempty(obj.hClust.inconsistentFields())
    % confirm recover with user
    if ~obj.hCfg.getOr('autoRecover', 0)
        dlgans = questdlg('Data found to be in an inconsistent state. Should I try to recover it?', 'Confirm Auto Recover', 'Yes');

        if ~strcmp(dlgans, 'Yes')
            success = 0;
            return;
        end
    end

    flag = obj.hClust.recover(1); % recover inconsistent data if needed
    successAppend = 'You should look through your data and ensure everything is correct, then save it.';
    failureAppend = 'You will probably experience problems curating your data.';
    msg = '';

    switch flag
        case 2
            msg = sprintf('Non-contiguous spike table found and corrected. %s', successAppend);

        case 1
            msg = sprintf('Inconsistent fields found and corrected. %s', successAppend);

        case 0
            msg = sprintf('Clustering data in an inconsistent state and automatic recovery failed. Please post an issue on the GitHub issue tracker. %s', failureAppend);
            success = 0;

        case -1
            msg = sprintf('Automatic recovery canceled by the user but the clustering data is still in an inconsistent state. %s', failureAppend);
            success = 0;
    end

    if ~isempty(msg)
        jrclust.utils.qMsgBox(msg, 1, 1);
    end
end
end %fun

