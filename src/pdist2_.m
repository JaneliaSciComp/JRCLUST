function Z = pdist2_(X, Y)
    % mr12 = pdist2_(mr1) % self distance
    % mr12 = pdist2_(mr1, mr2)
    % mr1: n1xd, mr2: n2xd, mr12: n1xn2

    % mr12 = sqrt(eucl2_dist_(mr1', mr2'));
    % 20% faster than pdist2 for 10000x10 x 100000x10 single
    % who cares how fast you are if you return complex numbers?
    if nargin < 2
        Y = X;
    end
    Z = pdist2(X, Y);

    return;

    if nargin == 2
        Z = sqrt(bsxfun(@plus, sum(Y'.^2), bsxfun(@minus, sum(X'.^2)', 2*(X*Y'))));
    else
        vr1 = sum(X'.^2);
        Z = sqrt(bsxfun(@plus, vr1, bsxfun(@minus, vr1', 2*(X*X'))));
    end
end
