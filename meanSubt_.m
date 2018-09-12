%--------------------------------------------------------------------------
% 8/6/17 JJJ: Tested and documented
function trWav1 = meanSubt_(trWav1, iDimm, hFunc)
    % subtract mean for mr or tr
    if nargin<2, iDimm = 1; end
    if nargin<3, hFunc = @mean; end
    if ~isa_(trWav1, 'single') && ~isa_(trWav1, 'double')
        trWav1 = single(trWav1);
    end
    trWav1_dim = size(trWav1);
    if numel(trWav1_dim)>2, trWav1 = reshape(trWav1, trWav1_dim(1), []); end
    trWav1 = bsxfun(@minus, trWav1, hFunc(trWav1, iDimm));
    if numel(trWav1_dim)>2, trWav1 = reshape(trWav1, trWav1_dim); end
end %func
