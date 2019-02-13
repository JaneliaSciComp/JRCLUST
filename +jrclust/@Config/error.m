function error(obj, emsg, varargin)
    %ERROR Raise an error
    obj.isError = 1;
    if obj.batchMode
        error(emsg);
    else
        uiwait(errordlg(emsg, varargin{:}, 'modal'));
    end
end