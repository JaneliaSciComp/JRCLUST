%--------------------------------------------------------------------------
function corr12 = cluWav_dist_(mrWav_clu1, mrWav_clu2)
    nPc = 3;
    viDelay = 0;

    [~, mrPv1, vrL1] = pca(mrWav_clu1, 'NumComponents', nPc);
    [~, mrPv2, vrL2] = pca(mrWav_clu2, 'NumComponents', nPc);
    % viDelay = -4:4; %account for time delay

    if viDelay==0
        corr12 = mean(abs(diag(corr_(mrPv1, mrPv2))));
    else
        vrDist12 = zeros(size(viDelay), 'single');
        for iDelay=1:numel(viDelay)
            mrPv2a = shift_mr_(mrPv2, viDelay(iDelay));
            vrDist12(iDelay) = mean(abs(diag(corr_(mrPv1, mrPv2a))));
        end
        corr12 = max(vrDist12);
    end
    % dist12 = abs(corr_(mrPv1, mrPv2));
    % [~,~,vrL12] = pca([mrWav_clu1, mrWav_clu2]);
    % nPc = 1;
    % dist12 = sum(vrL12(1:nPc)) / sum(vrL12);
    %
    % [~,~,vrL1] = pca(mrWav_clu1);
    % dist1 = sum(vrL1(1:nPc)) / sum(vrL1);
    % % dist1 = vrL(1) / sum(vrL);
    % %
    % [~,~,vrL2] = pca(mrWav_clu2);
    % dist2 = sum(vrL2(1:nPc)) / sum(vrL2);
    % % dist2 = vrL(1) / sum(vrL);
    % %
    % dist12 = dist12 / ((dist1 + dist2)/2);

    % mrCov12 = [mrWav_clu1, mrWav_clu2];
    % mrCov12 = mrCov12 * mrCov12';
    % [~,vrD] = eig([mrWav_clu1, mrWav_clu2]);
    % dist12 = vrD(end) / sum(vrD);

    % nPc = 3;
    % mrPc12 = eigvec_(mrCov12 * mrCov12', nPc);
    % mrPc1 = eigvec_(mrWav_clu1 * mrWav_clu1', nPc);
    % mrPc2 = eigvec_(mrWav_clu2 * mrWav_clu2', nPc);
    % dist12 = mean(abs(mean(mrPc2 .* mrPc1)));
end %func
