function depWarn(cmd, infoStr)
    %DEPWARN warn user of deprecated functionality
    % optionally, print a helpful message

    fprintf(2, '`%s` has been deprecated\n', cmd);
    if nargin > 1 && ~isempty(infoStr)
        fprintf(2, '%s\n', infoStr);
    end
end