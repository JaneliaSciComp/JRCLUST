function vals = rowColSelect(mat, rows, cols)
    %ROWCOLSELECT Select elements of a matrix by their row and column indices
    %   e.g. rowColSelect(eye(3), 1:3, 1:3) returns [1; 1; 1]
    if isempty(mat)
        vals = [];
        return;
    end

    if isempty(rows)
        rows = 1:size(mat,1);
    end

    if isempty(cols)
        cols = 1:size(mat, 2);
    end

    if numel(rows) ~= numel(cols)
        vals = [];
        return;
    end

    vals = mat(sub2ind(size(mat), rows(:), cols(:)));
end
