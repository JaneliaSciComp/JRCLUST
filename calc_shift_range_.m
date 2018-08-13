%--------------------------------------------------------------------------
% 10/27/17 JJJ: distance-based neighboring unit selection
function [cvi1, cvi2] = calc_shift_range_(P)

    % compute center range for the raw spike waveform
    spkLim_raw = get_(P, 'spkLim_raw');
    nSamples_raw = diff(spkLim_raw) + 1;
    spkLim_factor_merge = get_set_(P, 'spkLim_factor_merge', 1);
    spkLim_merge = round(P.spkLim * spkLim_factor_merge);
    viRange = (spkLim_merge(1) - spkLim_raw(1) + 1):(nSamples_raw - spkLim_raw(2) + spkLim_merge(2));

    % compute shift
    nShift = ceil(P.spkRefrac_ms / 1000 * P.sampleRateHz); % +/-n number of samples to compare time shift
    [cvi1, cvi2] = deal(cell(nShift*2+1, 1));
    viShift = -nShift:nShift;
    for iShift_ = 1:numel(viShift)
        iShift = viShift(iShift_);
        iShift1 = -round(iShift/2);
        iShift2 = iShift + iShift1;
        viRange1 = viRange + iShift1;
        viRange2 = viRange + iShift2;
        vl12 = (viRange1>=1 & viRange1<=nSamples_raw) & (viRange2>=1 & viRange2<=nSamples_raw);
        cvi1{iShift_} = viRange1(vl12);
        cvi2{iShift_} = viRange2(vl12);
    end
end %func
