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
        catch ME
            errMsg = sprintf('Could not set "%s": %s', propname, ME.message);
            flag = 0;
        end            
        try
            % transform val in some way
            if isfield(validData, 'postapply')
                hFun = eval(validData.postapply);
                val = hFun(val);
            end
        catch
            errMsg = sprintf('Could not set "%s" because postapply function (%s) failed.', propname, validData.postapply);
            flag = 0;
        end            
        try
            % check additional constraints
            if isfield(validData, 'postassert')
                hFun = eval(validData.postassert);
                assert(hFun(val));
            end
        catch
            errMsg = sprintf('Could not set "%s" because the following assertion (%s) failed.', propname, validData.postassert);
            flag = 0;
        end              
        try
            if isfield(validData, 'values')
                assert(all(ismember(val, validData.values)));
            end
        catch
            if iscellstr(validData.values)
                errMsg = sprintf('Could not set "%s" because it is not one of the required values: ', propname);                
                for i=1:length(validData.values)
                    if i<length(validData.values)
                        errMsg = [errMsg sprintf('"%s", ',validData.values{i})];
                    else
                        errMsg = [errMsg sprintf('"%s"',validData.values{i})];
                    end
                end
            else
                errMsg = sprintf('Could not set "%s" because it is not one of the required values.', propname);                
               % if list of required values isn't all strings, I haven't
               % tried to implement printing them all out. (TODO?)
            end
            flag = 0;
        end
    end
end