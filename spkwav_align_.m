%--------------------------------------------------------------------------
% 10/18/17 JJJ: created
function spikeWaveforms1 = spkwav_align_(spikeWaveforms1, P)
    nInterp_spk = 2;
    % viTime0 = 1:size(spikeWaveforms,1);

    imin_int = (-P.spkLim(1))*nInterp_spk+1;
    [mnWav_spk1_int, vi_int] = interpft_(spikeWaveforms1(:,:,1), nInterp_spk);
    [~,viMin_int] = min(mnWav_spk1_int);
    viSpk_right = find(viMin_int == imin_int+1);
    viSpk_left = find(viMin_int == imin_int-1);
    if ~isempty(viSpk_right)
        spikeWaveforms1_right = interpft_(spikeWaveforms1(:,viSpk_right,:), nInterp_spk);
        spikeWaveforms1(1:end-1,viSpk_right,:) = spikeWaveforms1_right(2:nInterp_spk:end,:,:);
    end
    if ~isempty(viSpk_left)
        spikeWaveforms1_left = interpft_(spikeWaveforms1(:,viSpk_left,:), nInterp_spk);
        spikeWaveforms1(2:end,viSpk_left,:) = spikeWaveforms1_left(2:nInterp_spk:end,:,:); %todo for nInterp_spk~=2
    end
end % function
