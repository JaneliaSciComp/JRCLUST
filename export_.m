%--------------------------------------------------------------------------
function export_(varargin)
    % export_(): export S0 struct to the workspace
    % export_(var1, var2): export fields in S0 struct to the workspace
    nArgs = nargin();
    csVars = varargin;
    csVars = csVars(~isempty_(csVars));
    if isempty(csVars), csVars = {'S0'}; end
    S0 = get(0, 'UserData');
    for iArg = 1:numel(csVars)
        vcVar = csVars{iArg};
        if isempty(vcVar), continue; end
        if ~strcmpi(vcVar, 'S0')
            var = get_(S0, vcVar);
        else
            var = S0;
        end
        if isempty(var)
            fprintf(2, '''%s'' does not exist\n', vcVar);
        else
            assignin('base', vcVar, var);
            fprintf('assigned ''%s'' to workspace\n', vcVar);
        end
    end
end % function
