%--------------------------------------------------------------------------
function trCorr = xcorr_mr_(mrPv_clu, nShift)
    % vrDist12 = xcorr_mr_(mrWav1, mrWav2, nShift)
    % vrDist12 = xcorr_mr_(mrWav1, mrWav2, cvi1, cvi2)
    if nShift==0
        trCorr = corr_(mrPv_clu);
        return;
    end
    nT = size(mrPv_clu,1);
    [cvi1, cvi2] = shift_range_(nT, nShift);
    nClu = size(mrPv_clu,2);
    trCorr = zeros([nClu, nClu, numel(cvi1)], 'like', mrPv_clu);
    for iShift = 1:numel(cvi1)
        trCorr(:,:,iShift) = corr_(mrPv_clu(cvi1{iShift},:), mrPv_clu(cvi2{iShift},:));
    end
end %func
