function prmFile = getPrmFile(filename)
    %GETPRMFILE Infer the .prm file from current working directory.
    %   Detailed explanation goes here

    if nargin < 1 % no filename specified
        d = dir('*.prm');
        if numel(d) ~= 1 % no .prm file or many found (ambiguous)
            error('Please specify a .prm file.');
        else
            prmFile = d(1).name;
        end
    else
        d = dir(filename);
        if numel(d) == 0 % filename not found
            error(['Filename "', filename, '" not found']);
        else
            prmFile = d(1).name;
        end
    end

end % func

