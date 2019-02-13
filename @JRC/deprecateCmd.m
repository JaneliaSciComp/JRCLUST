function deprecateCmd(obj, oldCmd, iMsg, newCmd)
    %DEPRECATECMD Print a message regarding deprecation
    %   if newCmd is specified, sets obj.cmd to newCmd
    if nargin < 3
        iMsg = '';
    end
    if nargin < 4
        newCmd = '';
    end

    % set a default iMsg in case of an alias
    if isempty(iMsg) && ~isempty(newCmd)
        iMsg = sprintf('Please use ''%s'' in the future', newCmd);
    end

    jrclust.utils.depWarn(oldCmd, iMsg);

    % set obj.cmd here
    if ~isempty(newCmd)
        obj.cmd = newCmd;
    end
end