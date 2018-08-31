%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function answer = userDialog(varargin)
    % Display a question dialog box
    global fDebug_ui
    if fDebug_ui == 1
        answer = 'Yes';
    else
        answer = questdlg(varargin{:});
    end
end % function
