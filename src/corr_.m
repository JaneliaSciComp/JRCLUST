%--------------------------------------------------------------------------
function C = corr_(A, B, fMeanSubt)
    % mr = corr_(A, B)
    % mr = corr_(A) % n1 x n2 becomes n1 x
    % mr = corr_(vr1, vr2) % single coefficient
    if nargin<3, fMeanSubt = 1; end
    % https://stackoverflow.com/questions/9262933/what-is-a-fast-way-to-compute-column-by-column-correlation-in-matlab
    if fMeanSubt, A = bsxfun(@minus,A,mean(A)); end %% zero-mean
    A = bsxfun(@times,A,1./sqrt(sum(A.^2))); %% L2-normalization
    if nargin == 1
        C = A' * A;
    else
        if fMeanSubt, B = bsxfun(@minus,B,mean(B)); end %%% zero-mean
        B = bsxfun(@times,B,1./sqrt(sum(B.^2))); %% L2-normalization
        C = A' * B;
    end
end %func
