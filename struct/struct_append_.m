%--------------------------------------------------------------------------
function S_old = struct_append_(S_old, S_new, varargin)
    if isempty(S_new), return; end
    if nargin==2
        varargin = fieldnames(S_new);
    end
    for i=1:numel(varargin)
        try
            S_old.(varargin{i}) = S_new.(varargin{i});
        catch
            ;
        end
    end %for
end % function
