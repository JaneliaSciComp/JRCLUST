%--------------------------------------------------------------------------
function [filteredTraces, vnWav1_mean] = wav_car_(filteredTraces, P)
    % take common average referencing (CAR) on the filtered trace (filteredTraces)
    % fprintf('Common average referencing (CAR)\n\t'); t1=tic;
    fRepairSites = 0; % bad sites get repaired by averaging vertical neighbors
    vnWav1_mean = [];
    switch lower(P.vcCommonRef)
        case {'tmean', 'nmean'}
            trimLim = [.25, .75];
            maxSite_ref = (P.nSites_ref + P.nSites_excl_ref - 1)/2;
            miSite_ref = findNearSites_(P.mrSiteXY, maxSite_ref, P.viSiteZero, P.viShank_site);
            miSite_ref = miSite_ref(P.nSites_excl_ref+1:end, :); %excl three nearest sites
            viChan_keep = round(trimLim * size(miSite_ref,1));
            viChan_keep = (viChan_keep(1)+1):viChan_keep(2);
            mnWav1_pre = filteredTraces;
            if strcmpi(P.vcCommonRef, 'tmean')
                for iChan=1:size(filteredTraces,2)
                    mnWav2 = sort(mnWav1_pre(:, miSite_ref(:,iChan)), 2);
                    gvr_tmean = sum(mnWav2(:,viChan_keep), 2); %may go out of range
                    gvr_tmean = int16(single(gvr_tmean)/numel(viChan_keep));
                    filteredTraces(:,iChan) = mnWav1_pre(:,iChan) - gvr_tmean;
                    fprintf('.');
                end
            else
                for iChan=1:size(filteredTraces,2)
                    gvr_tmean = sum(mnWav1_pre(:, miSite_ref(:, iChan)), 2); %may go out of range
                    gvr_tmean = int16(single(gvr_tmean)/size(miSite_ref,1));
                    filteredTraces(:,iChan) = mnWav1_pre(:,iChan) - gvr_tmean;
                    fprintf('.');
                end
            end
        case 'mean'
            vnWav1_mean = mean_excl_(filteredTraces, P);
            filteredTraces = bsxfun(@minus, filteredTraces, vnWav1_mean);
        case 'median'
            vnWav1_median = median_excl_(filteredTraces, P);
            filteredTraces = bsxfun(@minus, filteredTraces, vnWav1_median);
        case 'whiten'
            filteredTraces = whiten_(filteredTraces, P);
    end

    filteredTraces(:, P.viSiteZero) = 0; % TW do not repair with fMeanSite_drift with was not used
    % % viSiteZero should be treated carefully. try to repair using nearest sites?
    % if get_(P, 'fMeanSite_drift')
    %     filteredTraces = meanSite_drift_(filteredTraces, P);
    % elseif fRepairSites
    %     filteredTraces = meanSite_drift_(filteredTraces, P, P.viSiteZero);
    % else
    %     filteredTraces(:, P.viSiteZero) = 0;
    % end
    % fprintf('\n\ttook %0.1fs.\n', toc(t1));
end % function
