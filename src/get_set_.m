function val = get_set_(hCfg, fieldName, defaultVal)
    % set a value if field does not exist (empty)
    if isempty(hCfg)
        val = defaultVal;
        return;
    end

    if isa(hCfg, 'jrclust.Config')
        if isprop(hCfg, fieldName)
            val = hCfg.(fieldName);
        else
            val = defaultVal;
        end
    elseif ~isstruct(hCfg)
        val = [];
        fprintf(2, 'get_set_: %s must be a struct\n', inputname(1));
        return;
    else
        val = get_(hCfg, fieldName);
    end

    if isempty(val)
        val = defaultVal;
    end
end %func
