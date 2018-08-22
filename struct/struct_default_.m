%--------------------------------------------------------------------------
function P = struct_default_(P, csName, def_val)
    % Set to def_val if empty or field does not exist
    % set the field(s) to default val

    if ischar(csName), csName = {csName}; end
    for iField = 1:numel(csName)
        vcName = csName{iField};
        if ~isfield(P, vcName)
            P.(vcName) = def_val;
        elseif isempty(P.(vcName))
            P.(vcName) = def_val;
        end
    end
end %func
