%--------------------------------------------------------------------------
% 10/13/17 JJJ: Created. Realign the spikes at the min
function [viSpk_shift, viShift] = spkwav_shift_(trWav, shift_max, P)
    % trWav: nT x nSpk x nChans
    imid = 1 - P.spkLim(1);
    if P.fDetectBipolar
        [~, viMin] = max(abs(trWav(:,:,1)));
    else
        [~, viMin] = min(trWav(:,:,1));
    end
    viSpk_shift = find(viMin ~= imid);
    viShift = imid - viMin(viSpk_shift);
    vlKeep = abs(viShift) <= shift_max;
    viSpk_shift = viSpk_shift(vlKeep);
    viShift = viShift(vlKeep);
    % for iShift = -shift_max:shift_max
    %     if iShift==0, continue; end
    %     viSpk1 = viSpk_shift(viShift == iShift);
    %     if iShift<0
    %         trWav(1:end+iShift,viSpk1,:) = trWav(1-iShift:end,viSpk1,:);
    %     else
    %         trWav(iShift+1:end,viSpk1,:) = trWav(1:end-iShift,viSpk1,:);
    %     end
    % end
end %func
