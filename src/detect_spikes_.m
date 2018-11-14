%--------------------------------------------------------------------------
function [viTime_spk, vrAmp_spk, viSite_spk] = detect_spikes_(mnWav3, vnThresh_site, vlKeep_ref, P)
    % fMerge_spk = 1;
    fMerge_spk = get_set_(P, 'fMerge_spk', 1);

    [n1, nSites, ~] = size(mnWav3);
    [cviSpk_site, cvrSpk_site] = deal(cell(nSites,1));

    fprintf('\tDetecting spikes from each channel.\n\t\t'); t1=tic;
    % parfor iSite = 1:nSites
    for iSite = 1:nSites
        % Find spikes
        [viSpk11, vrSpk11] = spikeDetectSingle_fast_(mnWav3(:,iSite), P, vnThresh_site(iSite));
        fprintf('.');

        % Reject global mean
        if isempty(vlKeep_ref)
            cviSpk_site{iSite} = viSpk11;
            cvrSpk_site{iSite} = vrSpk11;
        else
            [cviSpk_site{iSite}, cvrSpk_site{iSite}] = select_vr_(viSpk11, vrSpk11, find(vlKeep_ref(viSpk11)));
        end
    end
    vnThresh_site = jrclust.utils.tryGather(vnThresh_site);
    nSpks1 = sum(cellfun(@numel, cviSpk_site));
    fprintf('\n\t\tDetected %d spikes from %d sites; took %0.1fs.\n', nSpks1, nSites, toc(t1));

    % Group spiking events using vrWav_mean1. already sorted by time
    if fMerge_spk
        fprintf('\tMerging spikes...'); t2=tic;
        [viTime_spk, vrAmp_spk, viSite_spk] = spikeMerge_(cviSpk_site, cvrSpk_site, P);
        fprintf('\t%d spiking events found; took %0.1fs\n', numel(viSite_spk), toc(t2));
    else
        viTime_spk = cell2mat_(cviSpk_site);
        vrAmp_spk = cell2mat_(cvrSpk_site);
        viSite_spk = cell2vi_(cviSpk_site);
        %sort by time
        [viTime_spk, viSrt] = sort(viTime_spk, 'ascend');
        [vrAmp_spk, viSite_spk] = multifun_(@(x)x(viSrt), vrAmp_spk, viSite_spk);
    end
    vrAmp_spk = jrclust.utils.tryGather(vrAmp_spk);

    % Group all sites in the same shank
    if get_set_(P, 'fGroup_shank', 0)
        [viSite_spk] = group_shank_(viSite_spk, P); % change the site location to the shank center
    end
end %func
