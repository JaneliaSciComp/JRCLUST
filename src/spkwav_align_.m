%--------------------------------------------------------------------------
% 10/18/17 JJJ: created
function tnWav_spk1 = spkwav_align_(tnWav_spk1, P)
    nInterp_spk = 2;
    % viTime0 = 1:size(tnWav_spk,1);

    imin_int = (-P.spkLim(1))*nInterp_spk+1;
    [mnWav_spk1_int, vi_int] = jrclust.utils.interpWindows(tnWav_spk1(:,:,1), nInterp_spk);
    [~,viMin_int] = min(mnWav_spk1_int);
    viSpk_right = find(viMin_int == imin_int+1);
    viSpk_left = find(viMin_int == imin_int-1);
    if ~isempty(viSpk_right)
        tnWav_spk1_right = jrclust.utils.interpWindows(tnWav_spk1(:,viSpk_right,:), nInterp_spk);
        tnWav_spk1(1:end-1,viSpk_right,:) = tnWav_spk1_right(2:nInterp_spk:end,:,:);
    end
    if ~isempty(viSpk_left)
        tnWav_spk1_left = jrclust.utils.interpWindows(tnWav_spk1(:,viSpk_left,:), nInterp_spk);
        tnWav_spk1(2:end,viSpk_left,:) = tnWav_spk1_left(2:nInterp_spk:end,:,:); %todo for nInterp_spk~=2
    end
end %func
