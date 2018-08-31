%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function vc = struct2str_(S)
    if ~isstruct(S), vc = ''; return; end

    csVars = fieldnames(S);
    vc = '';
    for iVar = 1:numel(csVars)
        vc = sprintf('%s%s = %s;', vc, csVars{iVar}, field2str_(S.(csVars{iVar})));
        if iVar<numel(csVars), vc = sprintf('%s\n', vc); end
    end %for
end % function
