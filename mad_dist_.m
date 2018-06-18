%--------------------------------------------------------------------------
function d12 = mad_dist_(mrFet1, mrFet2)
    % distance between two clusters
    if ~ismatrix(mrFet1)
        mrFet1 = reshape(mrFet1, [], size(mrFet1,3));
    end
    if ~ismatrix(mrFet2)
        mrFet2 = reshape(mrFet2, [], size(mrFet2,3));
    end
    vrFet1_med = median(mrFet1, 2);
    vrFet2_med = median(mrFet2, 2);
    vrFet12_med = vrFet1_med - vrFet2_med;
    norm12 = sqrt(sum(vrFet12_med.^2));
    vrFet12_med1 = vrFet12_med / norm12;
    mad1 = median(abs(vrFet12_med1' * bsxfun(@minus, mrFet1, vrFet1_med)));
    mad2 = median(abs(vrFet12_med1' * bsxfun(@minus, mrFet2, vrFet2_med)));
    d12 = norm12 / sqrt(mad1.^2 + mad2.^2);
end %func
