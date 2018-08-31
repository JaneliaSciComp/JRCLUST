function paramFile = getParamFile(filename)
    %GETparamFile Infer the .prm file from current working directory.
    %   Detailed explanation goes here

    if nargin < 1 % no filename specified
        d = dir('*.prm');
        if isempty(d)
            error('Please specify a .prm file.');
        else
            names = cell(numel(d), 1);
            for i=1:numel(d)
                names{i} = d(i).name;
            end
            notFull = ~endsWith(names, '_full.prm');
            if sum(notFull) == 1
                paramFile = names{notFull};
            else
                error('Ambiguity in param file; please specify.');
            end
        end
    else
        d = dir(filename);
        if numel(d) == 0 % filename not found
            error(['Filename "', filename, '" not found']);
        else
            paramFile = d(1).name;
        end
    end

end % function

