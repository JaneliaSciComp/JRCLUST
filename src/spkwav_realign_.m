%--------------------------------------------------------------------------
% 10/13/17 JJJ: Created. Realign the spikes at the min
function [tnWav_spk1, viTime_spk1] = spkwav_realign_(tnWav_spk1, mnWav_spk, spkLim_wav, viTime_spk1, viSite1, P)
    % subtract car and temporal shift
    % tnWav_spk1: nSamples x nSpk x nSites_spk
    if ~strcmpi(get_set_(P, 'vcSpkRef', 'nmean'), 'nmean'), return; end

    % fprintf('\n\tRealigning spikes after LCAR (vcSpkRef=nmean)...'); t1=tic;
    dimm1 = size(tnWav_spk1);
    trWav_spk2 = spkwav_car_(single(tnWav_spk1), P); % apply car

    shift_max = get_set_(S0.P, 'spkRefrac', []);
    if isempty(shift_max)
        shift_max = round(get_set_(S0.P, 'spkRefrac_ms', 0.25)*P.sRateHz/1000);
    end    
    [viSpk_shift, viShift] = spkwav_shift_(trWav_spk2, shift_max, P);

    [viSpk_shift, viShift] = spkwav_shift_(trWav_spk2, 1, P);
    if isempty(viSpk_shift), return; end

    viTime_shift = viTime_spk1(viSpk_shift) - int32(viShift(:)); % spike time to shift
    viTime_spk1(viSpk_shift) = viTime_shift;
    tnWav_spk1(:,viSpk_shift,:) = mr2tr3_(mnWav_spk, spkLim_wav, viTime_shift, viSite1);
    % fprintf('\n\t\ttook %0.1fs\n', toc(t1));
end %func
