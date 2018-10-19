%--------------------------------------------------------------------------
% 12/28/17 JJJ: P can be empty (v3.2.1)
% 8/4/17 JJJ: selective struct merge
% 7/31/17 JJJ: documentation and testing
function P = struct_merge_(P, P1, csNames)
    % Merge second struct to first one
    % P = struct_merge_(P, P_append)
    % P = struct_merge_(P, P_append, var_list) : only update list of variable names
    if isempty(P), P=P1; return; end % P can be empty
    if isempty(P1), return; end
    if nargin<3, csNames = fieldnames(P1); end
    if ischar(csNames), csNames = {csNames}; end

    for iField = 1:numel(csNames)
        vcName_ = csNames{iField};
        if isfield(P1, vcName_), P.(vcName_) = P1.(vcName_); end
    end
end %func
