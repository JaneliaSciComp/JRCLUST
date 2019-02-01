function mat_ = neCell2mat(cell_)
    %NECELL2MAT Like cell2mat, but keeps only nonempty cells
    nonempty = cellfun(@(x) ~isempty(x), cell_);
    mat_ = cell2mat(cell_(nonempty));
end
