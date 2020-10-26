function mat_ = neCell2mat(cell_)
    %NECELL2MAT Like cell2mat, but keeps only nonempty cells
    nonempty = cellfun(@(x) ~isempty(x), cell_);
    cell_ = cellfun(@(c) c(:), cell_(nonempty), 'UniformOutput', 0);

    mat_ = cell2mat(cell_(:));
end
