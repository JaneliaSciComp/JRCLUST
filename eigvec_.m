%--------------------------------------------------------------------------
function mrCov2 = eigvec_(mrCov2, nPc)
    [mrCov2,~] = eig(mrCov2);
    mrCov2 = mrCov2(:, end-nPc+1:end);
    mrCov2 = bsxfun(@minus, mrCov2, mean(mrCov2));
    mrCov2 = bsxfun(@rdivide, mrCov2, std(mrCov2));
end %func
