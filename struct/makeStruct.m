%--------------------------------------------------------------------------
function S = makeStruct(varargin)
    % MAKESTRUCT create a struct from varargin.
    % All inputs must be variables, not, e.g., functions of variables.

    S = struct();
    for i = 1:nargin
        S = setfield(S, inputname(i), varargin{i});
    end
end % function
