function filler = makeFiller(vals, shape)
%MAKEFILLER Make empty data similar to `vals` with shape `shape`.
cls = class(vals);
switch cls
    case 'cell' % expecting a homogeneous cell array, e.g., of char
        filler = cell(shape);
        if ~isempty(vals)
            filler = cellfun(@(x) cast(x, 'like', vals{1}), filler, 'UniformOutput', 0);
        end

    case 'char'
        filler = char(shape);

    otherwise
        filler = zeros(shape, cls);
end
end % func