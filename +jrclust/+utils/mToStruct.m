function S = mToStruct(filename)
    %MTOSTRUCT read and evaluate a .m script, convert workspace to struct
    lines = jrclust.utils.readLines(filename);
    lines = stripComments(lines);
    
    if isempty(lines)
        S = [];
        return;
    end

    try
        eval(cell2mat(lines'));
        wspace = whos();
        varnames = setdiff({wspace.name}, {'lines', 'filename'});

        for j = 1:numel(varnames)
            eval(sprintf('a = %s;', varnames{j}));
            S.(varnames{j}) = a;
        end
    catch
        S = struct();
    end
end

%% local functions
function lines = stripComments(lines)
    lines = lines(cellfun(@(ln) ~isempty(ln), lines));
    lines = cellfun(@(ln) strtrim(ln), lines, 'UniformOutput', 0);
    lines = lines(cellfun(@(ln) ln(1) ~= '%', lines));

    % remove comments in the middle
    for i = 1:numel(lines)
        iLine = lines{i};

        cMarker = find(iLine == '%', 1, 'first');
        if ~isempty(cMarker)
            iLine = iLine(1:cMarker - 1);
        end

        iLine = strrep(iLine, '...', '');
        if ismember(strsplit(iLine), {'for', 'end', 'if'})
            lines{i} = [strtrim(iLine), ', ']; % add blank at the end
        else
            lines{i} = [strtrim(iLine), ' ']; % add blank at the end
        end
    end

    lines = lines(cellfun(@(ln) ~isempty(ln), lines));
end