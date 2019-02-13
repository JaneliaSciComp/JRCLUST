function warning(obj, wmsg, varargin)
    %WARNING Raise a warning
    if obj.batchMode
        warning(wmsg);
    else
        uiwait(warndlg(wmsg, varargin{:}, 'modal'));
    end
end