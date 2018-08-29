%--------------------------------------------------------------------------
% 10/16/17 JJJ: Created
function [tr, vi] = interpft_(tr, nInterp)
    if nInterp==1, return; end
    vi = 1:(1/nInterp):size(tr,1);
    tr = interpft(tr, size(tr,1)*nInterp);
    if ndims(tr)==3
        tr = tr(1:numel(vi),:,:);
    else
        tr = tr(1:numel(vi),:);
    end
end % function
