%--------------------------------------------------------------------------
% 10/11/17 JJJ: mean_align_spkwav
function mrWav_mean = mean_align_spkwav_(tnWav_spk, P)
    tnWav_spk = meanSubt_(tnWav_spk);

    for iRepeat = 1:2
        mrWav_spk = single(reshape(tnWav_spk(:,1,:), [], size(tnWav_spk,3))); % use max chan only for selection
        % vrWav_mean = zscore(median(mrWav_spk,2), 1);
        vrCorr = zscore(mean(mrWav_spk,2))' * zscore(mrWav_spk); % / numel(vrWav_mean);
        %     viSpk_mean = find(vrCorr > quantile(vrCorr, 1/4));
        tnWav_spk = tnWav_spk(:,:, vrCorr > quantile(vrCorr, 1/4));
    end

    mrWav_mean = mean(tnWav_spk, 3);
end %func
