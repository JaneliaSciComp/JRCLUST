%--------------------------------------------------------------------------
function csAns = inputdlg_(varargin)
    % return default answer
    global fDebug_ui
    % if getOr([], 'fDebug_ui', 0)
    if fDebug_ui==1
        if numel(varargin)==4
            csAns = varargin{4};
        else
            csAns = [];
        end
    else
        csAns = inputdlg(varargin{:});
    end
end
