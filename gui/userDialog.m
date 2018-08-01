%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function vcAns = userDialog(varargin)
    % Display a question dialog box
    global fDebug_ui
    if fDebug_ui == 1
        vcAns = 'Yes';
    else
        vcAns = questdlg(varargin{:});
    end
end % func
