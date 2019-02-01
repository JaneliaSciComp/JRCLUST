function depWarn(cmd, infoStr)
    %DEPWARN warn user of deprecated functionality
    % optionally, print a helpful message

    wmsg = sprintf('''%s'' has been deprecated', cmd);
    if nargin > 1 && ~isempty(infoStr)
        wmsg = sprintf('%s\n\n%s', wmsg, infoStr);
    end
    
    warndlg(wmsg, 'Deprecation warning');
end