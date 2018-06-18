%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function vcAns = questdlg_(varargin)
    % Display a question dialog box
    global fDebug_ui
    % if get_set_([], 'fDebug_ui', 0)
    if fDebug_ui == 1
        vcAns = 'Yes';
    else
        vcAns = questdlg(varargin{:});
    end
end %func
