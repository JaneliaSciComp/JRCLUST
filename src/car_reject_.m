%--------------------------------------------------------------------------
function [vlKeep_ref, vrMad_ref] = car_reject_(vrWav_mean1, P)
    blank_period_ms = get_set_(P, 'blank_period_ms', 5);
    blank_thresh = get_set_(P, 'blank_thresh', []);
    [vlKeep_ref, vrMad_ref] = deal([]);
    % if isempty(blank_thresh) || blank_thresh==0, return; end

    % tbin_ref = .01; %10 msec bin
    vrWav_mean1 = single(vrWav_mean1);
    nwin = round(P.sRateHz * blank_period_ms / 1000);
    if nwin > 0 %nwin <= 1
        if nargout < 2
            vlKeep_ref = thresh_mad_(abs(vrWav_mean1), blank_thresh);
        else
            [vlKeep_ref, vrMad_ref] = thresh_mad_(abs(vrWav_mean1), blank_thresh);
        end

        if nwin > 1
            over_thresh=find(~vlKeep_ref);
            for crossing=over_thresh'
                left = max(crossing-ceil(nwin/2), 1);
                right = min(crossing+ceil(nwin/2), numel(vlKeep_ref));
                vlKeep_ref(left:right)=0;
            end
        end
    else
        vrRef_bin = std(reshape_vr2mr_(vrWav_mean1, nwin), 1,1);
        if nargout < 2
            vlKeep_ref = thresh_mad_(vrRef_bin, blank_thresh);
        else
            [vlKeep_ref, vrMad_ref] = thresh_mad_(vrRef_bin, blank_thresh);
            vrMad_ref = expand_vr_(vrMad_ref, nwin, size(vrWav_mean1));
        end
        vlKeep_ref = expand_vr_(vlKeep_ref, nwin, size(vrWav_mean1));
    end
    % figure; plot(vrMad_ref); hold on; plot(find(~vlKeep_ref), vrMad_ref(~vlKeep_ref), 'r.')
end %func
