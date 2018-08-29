%--------------------------------------------------------------------------
% 17/9/13 JJJ: Behavior changed, if S==[], S0 is loaded
function val = getOr(S, fieldName, defaultValue)
    % load a value from a struct or return a default if not found

    if isempty(S)
        S = get(0, 'UserData');
    end
    if isempty(S)
        val = defaultValue;
        return;
    end
    if ~isstruct(S)
        val = [];
        fprintf(2, 'getOr: %s must be a struct\n', inputname(1));
        return;
    end

    val = get_(S, fieldName);
    if isempty(val)
        val = defaultValue;
    end
end % function
