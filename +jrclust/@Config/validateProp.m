function [flag, val, errMsg] = validateProp(obj, propname, val)
    %VALIDATEPROP Ensure a property is valid
    if isfield(obj.oldParamSet, propname) % map the old param name to the new one
        propname = obj.oldParamSet.(propname);
        obj.isV3Import = 1;
    end

    flag = 1;
    errMsg = '';

    % found in default params, do validation
    if isfield(obj.fullParams, propname)
        validData = obj.fullParams.(propname).validation;
        classes = validData.classes;
        attributes = validData.attributes;

        if isempty(val) || isempty(attributes)
            if ~any(cellfun(@(c) isa(val, c), classes))
                flag = 0;
            end

            return;
        end

        try
            validateattributes(val, classes, attributes);

            % transform val in some way
            if isfield(validData, 'postapply')
                hFun = eval(validData.postapply);
                val = hFun(val);
            end

            % check additional constraints
            if isfield(validData, 'postassert')
                hFun = eval(validData.postassert);
                assert(hFun(val));
            end

            if isfield(validData, 'values')
                assert(all(ismember(val, validData.values)));
            end
        catch ME
            errMsg = sprintf('Could not set %s: %s', propname, ME.message);
            flag = 0;
        end
    end
end