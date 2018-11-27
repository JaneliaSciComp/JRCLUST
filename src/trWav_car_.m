%--------------------------------------------------------------------------
% CAR subtraction using outer half of sites
function tr = trWav_car_(tr, P)
    vcSpkRef = get_set_(P, 'vcSpkRef', 'nmean');
    if ~strcmpi(P.vcSpkRef, 'nmean'), return; end
    tr = single(permute(tr, [1,3,2]));
    dimm_tr = size(tr);
    viSite_ref = ceil(size(tr,3)/2):size(tr,3);
    mrWav_ref = mean(tr(:,:,viSite_ref), 3);
    tr = jrclust.utils.meanSubtract(reshape(bsxfun(@minus, reshape(tr,[],dimm_tr(3)), mrWav_ref(:)), dimm_tr));
    tr = permute(tr, [1,3,2]);
end %func
