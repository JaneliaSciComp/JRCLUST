%--------------------------------------------------------------------------
% 11/1/17 JJJ: Created
function tr1 = zero_start_(tr1)
    % subtract the first
    dimm1 = size(tr1);
    if numel(dimm1) ~= 2, tr1 = reshape(tr1, dimm1(1), []); end

    tr1 = bsxfun(@minus, tr1, tr1(1,:));

    if numel(dimm1) ~= 2, tr1 = reshape(tr1, dimm1); end
end %func
