%--------------------------------------------------------------------------
% 10/13/17 JJJ: Created
function [tnWav_spk, viSpk] = spkwav_filter_(tnWav_spk, qThresh)
    % vrCorr_spk = spkwav_maxcor_(tnWav_spk);
    vrCorr_spk = spkwav_maxcor_(tnWav_spk(:,1,:));
    viSpk = find(vrCorr_spk > quantile(vrCorr_spk, qThresh));
    tnWav_spk = tnWav_spk(:,:,viSpk);
end %func
