function strval = field2str(val)
    %FIELD2STR Convert a value to a canonical string representation
    switch class(val)
        case {'int', 'int16', 'int32', 'uint16', 'uint32'}
            formatstr = '%d';

        case {'double', 'single', 'logical'}
            if numel(val) == 1 && mod(val(1), 1) == 0
                formatstr = '%d';
            else
                formatstr = '%g';
            end

        case 'char'
            strval = sprintf('''%s''', val);
            return;

        case 'cell'
            strval = '{';
            for i = 1:numel(val)
                strval = [strval, jrclust.utils.field2str(val{i})];
                if i < numel(val)
                    strval = [strval, ', '];
                end
            end
            strval = [strval, '}'];
            return;

        case 'struct'
            strval = 'struct(';
            strFields = fieldnames(val);
            for i = 1:numel(strFields) - 1
                fn = strFields{i};
                strval = [strval, sprintf('''%s'', ', fn), ...
                          jrclust.utils.field2str(val.(fn)), ', '];
            end
            fn = strFields{end};
            strval = [strval, sprintf('''%s'', ', fn), ...
                          jrclust.utils.field2str(val.(fn)), ')'];
            return;

        otherwise
            error('Unsupported format: %s\n', class(val));
    end

    if numel(val) == 1 % scalar
        strval = sprintf(formatstr, val);
    else % matrix or array
        strval = '[';
        for iRow = 1:size(val, 1)
            for iCol=1:size(val, 2)
                strval = [strval, jrclust.utils.field2str(val(iRow, iCol))];
                if iCol < size(val, 2)
                    strval = [strval, ', '];
                end
            end

            if iRow < size(val, 1)
                strval = [strval, '; '];
            end
        end

        strval = [strval, ']'];
    end
end
