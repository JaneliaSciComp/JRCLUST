%--------------------------------------------------------------------------
% 17/12/5 JJJ: Error info is saved
% Display error message and the error stack
function disperr_(vcMsg, hErr)
    % disperr_(vcMsg): error message for user
    % disperr_(vcMsg, hErr): hErr: MException class
    % disperr_(vcMsg, vcErr): vcErr: error string
    try
        dbstack('-completenames'); % display an error stack
        if nargin<1, vcMsg = ''; end
        if nargin<2, hErr = lasterror('reset');  end
        if ischar(hErr) % properly formatted error
            vcErr = hErr;
        else
            save_err_(hErr, vcMsg); % save hErr object?
            vcErr = hErr.message;
        end
    catch
        vcErr = '';
    end
    if nargin==0
        fprintf(2, '%s\n', vcErr);
    elseif ~isempty(vcErr)
        fprintf(2, '%s:\n\t%s\n', vcMsg, vcErr);
    else
        fprintf(2, '%s:\n', vcMsg);
    end
    try gpuDevice(1); disp('GPU device reset'); catch, end
end %func
