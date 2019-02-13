function setProp(obj, propname, val)
    %SETPROP Set a property
    if isfield(obj.oldParamSet, propname)
        propname = obj.oldParamSet.(propname);
        obj.isV3Import = 1;
    end

    if isfield(obj.fullParams, propname)
        if ~isprop(obj, propname)
            obj.addprop(propname);
        end

        obj.(propname) = val;
    elseif ismember(propname, {'singleRaw', 'multiRaw'}) % separate validation for these
        obj.(propname) = val;
    else
        obj.setCustomProp(propname, val);
    end
end