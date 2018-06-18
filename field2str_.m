%---------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function vcStr = field2str_(val)
    % convert a value to a strong

    switch class(val)
        case {'int', 'int16', 'int32', 'uint16', 'uint32'}
        vcFormat = '%d';
        case {'double', 'single'}
        if numel(val)==1 && mod(val(1),1)==0
            vcFormat = '%d';
        else
            vcFormat = '%g';
        end
        case 'char'
        vcStr = sprintf('''%s''', val);
        return;
        case 'cell'
        vcStr = '{';
        for i=1:numel(val)
            vcStr = [vcStr, field2str_(val{i})];
            if i<numel(val), vcStr = [vcStr, ', ']; end
        end
        vcStr = [vcStr, '}'];
        return;
        otherwise
        vcStr = '';
        fprintf(2, 'field2str_: unsupported format: %s\n', class(val));
        return;
    end

    if numel(val) == 1
        vcStr = sprintf(vcFormat, val);
    else % Handle a matrix or array
        vcStr = '[';
        for iRow=1:size(val,1)
            for iCol=1:size(val,2)
                vcStr = [vcStr, field2str_(val(iRow, iCol))];
                if iCol < size(val,2), vcStr = [vcStr, ', ']; end
            end
            if iRow<size(val,1), vcStr = [vcStr, '; ']; end
        end
        vcStr = [vcStr, ']'];
    end
end %func
