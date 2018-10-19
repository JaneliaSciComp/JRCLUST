%--------------------------------------------------------------------------
function disp_score_(vrSnr, vrFp, vrFn, vrAccuracy, vnSite, vnSpk, fVerbose)
    P = get0_('P');
    snr_thresh_gt = get_set_(P, 'snr_thresh_gt', 7);
    if nargin<6, fVerbose = 1; end
    fSort_snr = 1;
    % if nargin<4, snr_thresh = 7; end %use 10 as a default
    % disp_score_(vrSnr, vrFp, vrFn, snr_thresh)
    fprintf('SNR>%d Groundtruth Units\n', snr_thresh_gt);
    viClu_gt = find(vrSnr > snr_thresh_gt);
    if fSort_snr
        [~, viSrt] = sort(vrSnr(viClu_gt), 'ascend');
        viClu_gt = viClu_gt(viSrt);
    end
    [vrSnr, vrFp_pct, vrFn_pct, vrAccuracy_pct] = ...
    multifun_(@(x)x(viClu_gt), vrSnr, vrFp*100, vrFn*100, vrAccuracy*100);

    vrScore2 = 100-vrFp_pct/2-vrFn_pct/2;
    fprintf('\tFalse Positive (%%): '); disp_stats_(vrFp_pct);
    fprintf('\tFalse Negative (%%): '); disp_stats_(vrFn_pct);
    fprintf('\tAccuracy (%%): '); disp_stats_(vrAccuracy_pct);
    fprintf('\tscore (1-FP-FN) (%%): '); disp_stats_(100-vrFp_pct-vrFn_pct);
    fprintf('\tscore2 1-(FP+FN)/2 (%%): '); disp_stats_(vrScore2);

    if nargin<=4, return; end
    [vnSite, vnSpk] = multifun_(@(x)x(viClu_gt), vnSite, vnSpk);
    [vrSnr, vrFp_pct, vrFn_pct, vrAccuracy_pct, vnSite, vnSpk, viClu_gt, vrScore2] = ...
    multifun_(@(x)x(:), vrSnr, vrFp_pct, vrFn_pct, vrAccuracy_pct, vnSite, vnSpk, viClu_gt, vrScore2);
    [vrSnr, vrFp_pct, vrFn_pct, vrAccuracy_pct, vrScore2] = ...
    multifun_(@(x)round(x*10)/10, vrSnr, vrFp_pct, vrFn_pct, vrAccuracy_pct, vrScore2);

    if fVerbose
        disp(table(viClu_gt, vrSnr, vrScore2, vrFp_pct, vrFn_pct, vrAccuracy_pct, vnSite, vnSpk));
    end
end %func
