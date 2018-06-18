%--------------------------------------------------------------------------
% 10/11/17 JJJ: Created
function [vrCorr_spk, vrCorr_spk2] = spkwav_maxcor_(tnWav_spk, tnWav_spk2)
    % subselect spikes that matches the mean well
    % dimm1 = size(tnWav_spk);
    % tnWav_spk0 = tnWav_spk;
    if nargin<2, tnWav_spk2 = []; end
    tnWav_spk = meanSubt_(single(gpuArray_(tnWav_spk)));
    mrWav_spk = zscore(reshape(tnWav_spk, [], size(tnWav_spk,3)), 1); % use max chan only for selection
    if isempty(tnWav_spk2)
        mrCorr_spk = set_diag_(mrWav_spk' * mrWav_spk, nan);
    else
        tnWav_spk2 = meanSubt_(single(gpuArray_(tnWav_spk2)));
        mrWav_spk2 = zscore(reshape(tnWav_spk2, [], size(tnWav_spk2,3)), 1); % use max chan only for selection
        mrCorr_spk = mrWav_spk2' * mrWav_spk;
    end
    vrCorr_spk = gather_(max(mrCorr_spk) / size(mrWav_spk,1));
    if nargout>1
        vrCorr_spk2 = gather_(max(mrCorr_spk,[],2) / size(mrWav_spk2,1));
    end

    % vrWav_mean = zscore(median(mrWav_spk,2), 1);
    %     vrCorr = zscore(mean(mrWav_spk,2))' * zscore(mrWav_spk); % / numel(vrWav_mean);
    % %     viSpk_mean = find(vrCorr > quantile(vrCorr, 1/4));
    %     tnWav_spk = tnWav_spk(:,:, vrCorr > quantile(vrCorr, 1/4));
    % end
    %
    % mrWav_mean = mean(tnWav_spk, 3);

    % [a,b,c] = pca(single(reshape(tnWav_spk(:,1,:), dimm1(1), [])), 'NumComponents', 1);
    % rankorder_(a);
    %
    % tnWav_spk0 = tnWav_spk;
    %
    % mrWav_spk = single(permute(tnWav_spk(:,1,:), [1,3,2]));
    % vrWav0 = median(mrWav_spk, 2);
    % nShift = numel(cviShift1);
    % nSpk = size(tnWav_spk,3);
    % mrCorr = zeros(nSpk, nShift);
    % % viShift = (1:nShift) - round(nShift/2);
    % for iShift = 1:nShift
    %     vrWav_ = zscore(vrWav0(cviShift1{iShift}), 1);
    %     mrWav_ = mrWav_spk(cviShift2{iShift}, :);
    %     mrCorr(:,iShift) = vrWav_' * mrWav_;
    % end
    % [~, viShift] = max(mrCorr, [], 2);
    % viShift = round(nShift/2) - viShift;
    % for iSpk = 1:nSpk
    %     iShift_ = viShift(iSpk);
    %     if iShift_ == 0, continue; end
    %     tnWav_spk(:,:,iSpk) = shift_mr_(tnWav_spk(:,:,iSpk), iShift_);
    % end
    % mrWav_mean = single(mean(tnWav_spk, 3));
end %func
