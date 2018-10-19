%--------------------------------------------------------------------------
% 17/12/11 JJJ: Created. Apply spatial whitening
function [mnWav2, mrWhiten] = whiten_(mnWav1, P)
    nLoads_gpu = get_set_(P, 'nLoads_gpu', 8);
    nSamples_max = round(size(mnWav1,1) / nLoads_gpu);
    fprintf('Whitening\n\t'); t1 = tic;
    [mr_sub, vi_sub] = subsample_mr_(mnWav1, nSamples_max, 1);
    viSites = setdiff(1:size(mnWav1,2), P.viSiteZero);
    if ~isempty(P.viSiteZero), mr_sub = mr_sub(:,viSites); end
    mr_sub = single(mr_sub);

    mrXXT = mr_sub' * mr_sub;
    [U,D] = eig(mrXXT + eps('single'));
    Sinv = diag(1./sqrt(diag(D)));
    scale = mean(sqrt(diag(mrXXT)));
    mrWhiten = (U * Sinv * U') * scale;

    % mr_sub1 = mr_sub * (U * Sinv * U');
    % mrXXT1 = mr_sub1' * mr_sub1;
    % figure; imagesc(mrXXT1 - eye(size(mrXXT1)))

    % apply whitening matrix
    mnWav2 = zeros(size(mnWav1), 'like', mnWav1);
    if ~isempty(P.viSiteZero), mnWav1 = mnWav1(:,viSites); end
    mnWav1 = single(mnWav1);
    for iSite1 = 1:numel(viSites)
        iSite = viSites(iSite1);
        mnWav2(:,iSite) = int16(mnWav1 * mrWhiten(:,iSite1));
        fprintf('.');
    end %for
    fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func
