%--------------------------------------------------------------------------
% 10/18/17 JJJ: speed optimization. Todo: use shift matrix?  (https://en.wikipedia.org/wiki/Shift_matrix)
% 10/8/17 JJJ: find correlation
function vrDist12 = xcorr2_mr_(mrWav1, mrWav2, arg1, arg2)
    % vrDist12 = xcorr_mr_(mrWav1, mrWav2, nShift)
    % vrDist12 = xcorr_mr_(mrWav1, mrWav2, cvi1, cvi2)
    % fMeanSubt_post = 1;
    % fSquared = 0;
    if nargin == 3
        nShift = arg1;
        nT = size(mrWav1, 1);
        [cvi1, cvi2] = shift_range_(nT, nShift);
    else
        cvi1 = arg1;
        cvi2 = arg2;
    end
    % if fSquared
    %     mrWav1 = mrWav1 .^ 2;
    %     mrWav2 = mrWav2 .^ 2;
    % end
    % vrDist12 = jrclust.utils.tryGpuArray(zeros(size(cvi1)), isGpu_(mrWav1));
    vrDist12 = zeros(size(cvi1));
    for iDist = 1:numel(vrDist12)
        vr1 = mrWav1(cvi1{iDist},:);
        vr2 = mrWav2(cvi2{iDist},:);
        %     vrDist12(iDist) = corr_(vr1(:), vr2(:), 1);
        vr1 = vr1(:);
        vr2 = vr2(:);
        %     n = numel(vr1);
        vr1 = vr1 - sum(vr1)/numel(vr1);
        vr2 = vr2 - sum(vr2)/numel(vr2);
        vrDist12(iDist) = vr1'*vr2 / sqrt(sum(vr1.^2)*sum(vr2.^2));
        %     vrDist12(iDist) = vr1'*vr2 / sqrt(vr1'*vr1*vr2'*vr2);
    end
end %func
