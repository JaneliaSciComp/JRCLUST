%--------------------------------------------------------------------------
function d = eucl2_dist_(X, Y)
    % a: m x d1; b: m x d2
    % aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b;
    % d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
    % X = [mrFet1_, mrFet2_];
    d = bsxfun(@plus, sum(Y.^2), bsxfun(@minus, sum(X.^2)', 2*X'*Y));
end % function
