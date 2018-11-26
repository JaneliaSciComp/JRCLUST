function z = zscore(X, w, dim)
    %ZSCORE Compute z-scores ((X - mean(X))/std(X)) for observations in X
    if isempty(X)
        z = [];
        return;
    end

    if nargin < 2
        w = 0;
    end
    if nargin < 3 % try to infer which dimension to work along
        dim = find(size(X) ~= 1, 1);
        if isempty(dim)
            dim = 1;
        end
    end

    mu = mean(X, dim);
    sigma = std(X, w, dim);
    sigma(sigma == 0) = 1; % no dispersion in sample, don't divide by zero

    z = bsxfun(@minus, X, mu);
    z = bsxfun(@rdivide, z, sigma);
end
