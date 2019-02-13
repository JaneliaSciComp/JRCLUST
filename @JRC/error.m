function error(obj, emsg)
    %ERROR Raise an error
    obj.isError = 1;
    obj.errMsg = emsg;
    error(emsg);
end