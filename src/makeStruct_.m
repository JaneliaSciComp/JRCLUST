%--------------------------------------------------------------------------
function S = makeStruct_(varargin)
    %MAKESTRUCT all the inputs must be a variable.
    %don't pass function of variables. ie: abs(X)
    %instead create a var AbsX an dpass that name
    S = struct();
    for i=1:nargin, S.(inputname(i)) =  varargin{i}; end
end %func
