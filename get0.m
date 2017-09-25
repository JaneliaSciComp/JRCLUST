function varargout = get0(varargin)
S0 = get(0, 'UserData'); 
if nargin==0
    varargout{1} = S0; 
    if nargout==0
        assignWorkspace(S0);
    end
    return; 
end

for i=1:nargin
    try
        
        if nargout==0
            eval(sprintf('%s = S0.%s;', varargin{i}, varargin{i}));
            eval(sprintf('assignWorkspace(%s);', varargin{i}));
%             eval(sprintf('clear %s;', varargin{i}));
%             fprintf('Assigned %s to workspace\n', varargin{i});
        else
            varargout{i} = S0.(varargin{i});
        end
    catch
        varargout{i} = [];
    end
end
end %func