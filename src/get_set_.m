%--------------------------------------------------------------------------
% 17/9/13 JJJ: Behavior changed, if S==[], S0 is loaded
function val = get_set_(S, vcName, def_val)
    % set a value if field does not exist (empty)

    if isempty(S), S = get(0, 'UserData'); end
    if isempty(S), val = def_val; return; end
    if isa(S, 'jrclust.Config')
        try
            val = S.(vcName);
        catch
            val = [];
        end
    elseif ~isstruct(S)
        val = [];
        fprintf(2, 'get_set_: %s must be a struct\n', inputname(1));
        return;
    else
        val = get_(S, vcName);
    end
    if isempty(val)
        val = def_val;
    end
end %func
