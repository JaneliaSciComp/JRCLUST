function d = eucl2_dist_(X, Y) % THIS IS JUST pdist2(X', Y')
    Z = bsxfun(@minus, sum(X.^2)', 2*X'*Y);
    d = bsxfun(@plus, sum(Y.^2), Z);
end
