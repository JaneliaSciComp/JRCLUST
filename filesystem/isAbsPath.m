function flag = isAbsPath(path)
    % ISABSPATH check whether a path given exists and is an absolute path
    if ~ischar(path)
        flag = 0;
    else
        if startsWith(lower(getenv('OS')), 'windows')
            flag = ~isempty(regexpi(path, '^[a-z]:')) && ~isempty(dir(path));
        else
            flag = startsWith(path, '/') && ~isempty(dir(path));
        end
    end
end

