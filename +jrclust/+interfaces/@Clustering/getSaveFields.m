function fieldNames = getSaveFields(obj)
%GETSAVEFIELDS Return a list of fields to save to or load from res file.
%   Save only non-dependent/non-transient fields, also ignore large fields
%   loaded from binary files.
m = metaclass(obj);
fieldNames = {m.PropertyList.Name};

% remove dependent, hidden, or transient fields
isDependent = [m.PropertyList.Dependent];
isHidden = [m.PropertyList.Hidden];
isTransient = [m.PropertyList.Transient];

fieldNames = fieldNames(~(isDependent | isTransient | isHidden));

end
