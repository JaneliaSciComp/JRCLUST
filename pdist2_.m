%--------------------------------------------------------------------------
function mr12 = pdist2_(mr1, mr2)
    % mr12 = pdist2_(mr1) % self distance
    % mr12 = pdist2_(mr1, mr2)
    % mr1: n1xd, mr2: n2xd, mr12: n1xn2

    % mr12 = sqrt(eucl2_dist_(mr1', mr2'));
    % 20% faster than pdist2 for 10000x10 x 100000x10 single
    if nargin==2
        mr12 = sqrt(bsxfun(@plus, sum(mr2'.^2), bsxfun(@minus, sum(mr1'.^2)', 2*mr1*mr2')));
    else
        vr1 = sum(mr1'.^2);
        mr12 = sqrt(bsxfun(@plus, vr1, bsxfun(@minus, vr1', 2*mr1*mr1')));
    end
end %func
