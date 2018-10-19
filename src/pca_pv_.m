%--------------------------------------------------------------------------
function [vrPv1, vrPv2] = pca_pv_(mr1)
    MAX_SAMPLE = 10000; %for pca
    mrCov = subsample_mr_(mr1, MAX_SAMPLE, 2);
    mrCov = meanSubt_(single(mrCov));
    mrCov = mrCov * mrCov';
    [mrPv1, vrD1] = eig(mrCov);
    vrPv1 = mrPv1(:,end);
    vrPv2 = mrPv1(:,end-1);
    % iMid = 1-P.spkLim(1);
    % vrPv1 = ifeq_(vrPv1(iMid)>0, -vrPv1, vrPv1);
    % vrPv2 = ifeq_(vrPv2(iMid)>0, -vrPv2, vrPv2);
end %func
