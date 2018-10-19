%--------------------------------------------------------------------------
function plot_activity_(P) % single column only
    % plot activity as a function of depth and time
    global fDebug_ui

    tbin = 10; %activity every 10 sec
    % plot activity as a function of time
    % vcFile_evt = subsFileExt(P.vcFile_prm, '_evt.mat');
    S0 = load_cached_(P, 0); %do not load waveforms
    nSites = numel(P.viSite2Chan);
    % tdur = max(cell2mat_(cellfun(@(x)double(max(x)), Sevt.cviSpk_site, 'UniformOutput', 0))) / P.sRateHz;
    tdur = double(max(S0.viTime_spk)) / P.sRateHz; % in sec
    nTime = ceil(tdur / tbin);

    mrAmp90 = zeros(nTime, nSites);
    lim0 = [1, tbin * P.sRateHz];
    for iSite=1:nSites
        viSpk1 = find(S0.viSite_spk == iSite);
        vrAmp_spk1 = S0.vrAmp_spk(viSpk1); % % S0.cvrSpk_site{iSite};  %spike amplitude
        if isempty(vrAmp_spk1), continue; end
        viTime_spk1 = S0.viTime_spk(viSpk1);
        for iTime=1:nTime
            lim1 = lim0 + (iTime-1) * lim0(2);
            vrAmp_spk11 = vrAmp_spk1(viTime_spk1 >= lim1(1) & viTime_spk1 <= lim1(2));
            if isempty(vrAmp_spk11),  continue; end
            mrAmp90(iTime, iSite) = quantile(abs(vrAmp_spk11), .9);
        end
    end %for
    mrAmp90=mrAmp90';

    vlSite_left = P.mrSiteXY(:,1) == 0;
    vrSiteY = P.mrSiteXY(:,2);
    hFig = create_figure_('FigActivity', [0 0 .5 1], P.vcFile_prm, 1, 1);
    subplot 121; imagesc(mrAmp90(vlSite_left, :), 'XData', (1:nTime) * tbin, 'YData', vrSiteY(vlSite_left)); axis xy; xlabel('Time'); ylabel('Sites'); title('Left edge sites');
    subplot 122; imagesc(mrAmp90(~vlSite_left, :), 'XData', (1:nTime) * tbin, 'YData', vrSiteY(~vlSite_left)); axis xy; xlabel('Time'); ylabel('Sites'); title('Right edge sites');

    [~, iSite_center] = max(mean(mrAmp90,2));
    viSiteA = iSite_center + [-2:2]; %look neighbors
    mrAmp90a = mrAmp90(viSiteA, :);
    vrCentroid = bsxfun(@rdivide, sum(bsxfun(@times, mrAmp90a.^2, vrSiteY(viSiteA))), sum(mrAmp90a.^2));
    hold on; plot((1:nTime) * tbin, vrCentroid, 'r');

    % if get_set_([], 'fDebug_ui', 0), close(hFig); end
    if fDebug_ui==1, close(hFig); end
end %func
