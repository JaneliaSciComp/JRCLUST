%--------------------------------------------------------------------------
% 9/26/17 JJJ: Output message is added
% 8/2/17 JJJ: Test and documentation
function vcMsg = assignWorkspace_(varargin)
    % Assign variables to the Workspace
    vcMsg = {};
    for i=1:numel(varargin)
        if ~isempty(varargin{i})
            assignin('base', inputname(i), varargin{i});
            vcMsg{end+1} = sprintf('assigned ''%s'' to workspace\n', inputname(i));
        end
    end
    vcMsg = cell2mat(vcMsg);
    if nargout==0, fprintf(vcMsg); end
end %func
