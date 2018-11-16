% function [vlKeep_ref, vrMad_ref] = carReject(vrWav_mean1, P)
function [vlKeep_ref, vrMad_ref] = carReject(channelMeans, blankPer_ms, blankThresh, sampleRate)
    %CARREJECT

    vrMad_ref = [];

    channelMeans = single(channelMeans);
    nwin = round(sampleRate * blankPer_ms / 1000);
    if nwin > 0 % nwin <= 1
        if nargout < 2
            vlKeep_ref = threshMAD(abs(channelMeans), blankThresh);
        else
            [vlKeep_ref, vrMad_ref] = threshMAD(abs(channelMeans), blankThresh);
        end

        if nwin > 1
            over_thresh = find(~vlKeep_ref);
            for crossing = over_thresh'
                left = max(crossing-ceil(nwin/2), 1);
                right = min(crossing+ceil(nwin/2), numel(vlKeep_ref));
                vlKeep_ref(left:right) = 0;
            end
        end
    else
        vrRef_bin = std(reshape_vr2mr_(channelMeans, nwin), 1,1);
        if nargout < 2
            vlKeep_ref = threshMAD(vrRef_bin, blankThresh);
        else
            [vlKeep_ref, vrMad_ref] = threshMAD(vrRef_bin, blankThresh);
            vrMad_ref = expand_vr_(vrMad_ref, nwin, size(channelMeans));
        end
        vlKeep_ref = expand_vr_(vlKeep_ref, nwin, size(channelMeans));
    end
end

%% LOCAL FUNCTIONS
function [vl, vr] = threshMAD(vr, thresh)
    %THRESHMAD single sided, no absolute value
    nSubs = 300000;
    offset = median(subsample_vr_(vr, nSubs));
    vr = vr - offset; % center the mean
    factor = median(abs(subsample_vr_(vr, nSubs)));

    if isempty(thresh) || thresh == 0
        vl = true(size(vr));
    else
        vl = vr < factor * thresh;
    end

    if nargout >= 2
        vr = vr / factor; % MAD unit
    end
end

function mr = reshape_vr2mr_(vr, nwin)
    nbins = ceil(numel(vr)/nwin);
    vr(nbins*nwin) = 0; % expand size
    mr = reshape(vr(1:nbins*nwin), nwin, nbins);
end
