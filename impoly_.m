%--------------------------------------------------------------------------
function hPoly = impoly_(varargin)
    global fDebug_ui
    % if get_set_([], 'fDebug_ui', 0)
    if fDebug_ui==1
        hPoly = []; %skip the test if debugging
    else
        hPoly = impoly(varargin{:});
    end
end
