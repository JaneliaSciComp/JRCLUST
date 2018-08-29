%--------------------------------------------------------------------------
function z = zscore_(x, flag, dim)
    if isempty(x), z=[]; return; end
    if nargin < 2, flag = 0; end
    if nargin < 3
        % Figure out which dimension to work along.
        dim = find(size(x) ~= 1, 1);
        if isempty(dim), dim = 1; end
    end

    % Compute X's mean and sd, and standardize it
    mu = mean(x,dim);
    sigma = std(x,flag,dim);
    sigma0 = sigma;
    sigma0(sigma0==0) = 1;
    z = bsxfun(@minus,x, mu);
    z = bsxfun(@rdivide, z, sigma0);
end % function
