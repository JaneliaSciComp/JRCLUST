%--------------------------------------------------------------------------
function vc = set_(vc, varargin)
    % Set handle to certain values
    % set_(S, name1, val1, name2, val2)

    if isempty(vc), return; end
    if isstruct(vc)
        for i=1:2:numel(varargin)
            vc.(varargin{i}) = varargin{i+1};
        end
        return;
    end
    if iscell(vc)
        for i=1:numel(vc)
            try
                set(vc{i}, varargin{:});
            catch
            end
        end
    elseif numel(vc)>1
        for i=1:numel(vc)
            try
                set(vc(i), varargin{:});
            catch
            end
        end
    else
        try
            set(vc, varargin{:});
        catch
        end
    end
end % function
