function mat = setDiag(mat, vec)
    %SETDIAG Set diagonal values of a matrix in a robust way
    n = min(min(size(mat)), numel(vec));
    mat(sub2ind(size(mat), 1:n, 1:n)) = vec(1:n);
end
