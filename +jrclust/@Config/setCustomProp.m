function setCustomProp(obj, propname, val)
    %SETCUSTOMPROP Set a property not included in the defaults
    if ismember(propname, obj.deprecatedParams.unsupported)
        return;
    end

    % ignore property if it's Dependent
    propData = ?jrclust.Config;
    propNames = {propData.PropertyList.Name};
    dependentProps = propNames([propData.PropertyList.Dependent]);
    if ismember(propname, dependentProps)
        return;
    end

    if ~isprop(obj, propname)
        addprop(obj, propname);
        if ~ismember(propname, obj.customParams)
            obj.customParams{end+1} = propname;
        end
    end

    obj.(propname) = val;
end