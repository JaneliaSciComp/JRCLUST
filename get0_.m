%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function varargout = get0_(varargin)
    % returns get(0, 'UserData') to the workspace
    % [S0, P] = get0_();
    % [S0, P, S_clu] = get0_();
    % [var1, var2] = get0_('var1', 'var2'); %sets [] if doesn't exist
    S0 = get(0, 'UserData');
    if ~isfield(S0, 'S_clu'), S0.S_clu = []; end
    if nargin==0
        varargout{1} = S0;
        if nargout==0, assignWorkspace_(S0); return; end
        if nargout>=1, varargout{1} = S0; end
        if nargout>=2, varargout{2} = S0.P; end
        if nargout>=3, varargout{3} = S0.S_clu; end
        return;
    end
    for i=1:nargin
        try
            eval(sprintf('%s = S0.%s;', varargin{i}, varargin{i}));
            varargout{i} = S0.(varargin{i});
        catch
            varargout{i} = [];
        end
    end
end % function
