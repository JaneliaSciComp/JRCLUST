%--------------------------------------------------------------------------
function hRect = imrect_(varargin)
    global fDebug_ui

    hRect = []; %skip the test if debugging
    % if get_set_([], 'fDebug_ui', 0) && nargin < 2
    if fDebug_ui==1 && nargin < 2
        return;
    else
        try
            hRect = imrect(varargin{:});
        catch
            fprintf(2, 'Install image processing toolbox\n');
        end
    end
end
