%--------------------------------------------------------------------------
function plot_drift_(P)

    iShank_show = getOr(P, 'iShank_show', 1);
    vcMode_drift = getOr(P, 'vcMode_drift', 'z'); % {'tay', 'xy', 'xya', 'x', 'y'}
    vcMode_com = getOr(P, 'vcMode_com', 'fet'); % {'fet', 'filt', 'raw', 'std'}
    % Compute spike position and plot drift profile
    S0 = load_cached_(P); % load cached data or from file if exists
    [S_clu, spikeSecondarySites, vrTime_spk, spikeSites] = get0_('S_clu', 'spikeSecondarySites', 'spikeTimes', 'spikeSites');
    vrTime_spk = double(vrTime_spk) / P.sampleRateHz;
    nSites_spk = 1 + P.maxSite*2 - P.nSites_ref;
    miSites_spk = P.miSites(:,spikeSites);

    switch lower(vcMode_com)
        case 'filt'
        spikeWaveforms = getSpikeWaveforms(P, 0);
        %mrVp = squeeze_(single(min(spikeWaveforms))) .^ 2;
        %         spikeWaveforms1 = meanSubt_(single(spikeWaveforms),2);
        mrVp = single(squeeze_(max(spikeWaveforms) - min(spikeWaveforms))) .^ 2;
        case 'filtstd'
        mrVp = squeeze_(var(single(getSpikeWaveforms(P, 0))));
        case 'rawstd', mrVp = squeeze_(var(single(getSpikeWaveforms(P, 1))));
        case 'raw'
        spikeTraces = getSpikeWaveforms(P, 1);
        mrVp = single(squeeze_(max(spikeTraces) - min(spikeTraces))) .^ 2;
        case 'fet'
        spikeFeatures = getSpikeFeatures(P);
        mrVp = squeeze_(spikeFeatures(1:nSites_spk,1,:)) .^ 2;
        %         mrVp = abs(squeeze_(spikeFeatures(1:nSites_spk,1,:)));
        miSites_spk = miSites_spk(1:nSites_spk,:);
    end

    mrX_spk = reshape(P.mrSiteXY(miSites_spk,1), size(miSites_spk));
    mrY_spk = reshape(P.mrSiteXY(miSites_spk,2), size(miSites_spk));
    vrPosX_spk = sum(mrVp .* mrX_spk) ./ sum(mrVp);
    vrPosY_spk = sum(mrVp .* mrY_spk) ./ sum(mrVp);
    vrA_spk = sum(mrVp);
    assignWorkspace_(vrTime_spk, vrPosX_spk, vrPosY_spk, vrA_spk);

    vlSpk_shank = ismember(P.viShank_site(spikeSites), iShank_show); %show first shank only
    % vrAmp_spk = 1 ./ sqrt(single(abs(S0.vrAmp_spk)));
    % vrAmp_spk = 1 ./ sqrt(sum(mrVp));
    vrAmp_spk = sqrt(mean(mrVp) ./ std(mrVp)); %spatial icv
    hFig_drift = createFigure('', [0 0 .5 1], P.paramFile, 1, 1);
    % hFig_drift = gcf;
    figure(hFig_drift);
    ax = gca();
    hold on;
    nSpk_thresh_clu = median(S_clu.nSpikesPerCluster);
    snr_thresh_clu = getOr(P, 'snr_thresh_clu', quantile(S_clu.vrSnr_clu, .5));
    % posX_thresh = median(S_clu.clusterXPositions);
    % [vrTime_drift, vrDepth_drift] = drift_track_(S_clu, vrPosY_spk, P);

    if ~isempty(S_clu)
        posX_lim = quantile(S_clu.clusterXPositions, [.25, .75]);
        posY_lim = quantile(S_clu.clusterYPositions, [.1, .9]);
        for iClu = 1:S_clu.nClusters
            viSpk1 = S_clu.spikesByCluster{iClu};
            viSpk1 = viSpk1(vlSpk_shank(viSpk1));
            vrColor1 = rand(1,3);
            %         vrColor1 = 'r';
            switch vcMode_drift
                case 'tay'
                plot3(ax, vrTime_spk(viSpk1), vrAmp_spk(viSpk1), vrPosY_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5);
                case 'x'
                plot(ax, vrTime_spk(viSpk1), vrPosX_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5);
                case 'y'
                plot(ax, vrTime_spk(viSpk1), vrPosY_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5);
                case 'xy'
                plot(ax, vrPosX_spk(viSpk1), vrPosY_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5);
                case 'xya'
                if S_clu.vrSnr_clu(iClu) < snr_thresh_clu, continue; end
                plot3(ax, vrPosX_spk(viSpk1), vrPosY_spk(viSpk1), vrAmp_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5);
                case 'z'
                %                 if S_clu.nSpikesPerCluster(iClu) < nSpk_thresh_clu, continue; end
                if S_clu.vrSnr_clu(iClu) < snr_thresh_clu, continue; end
                %                 if vrA_spk(S_clu.spikesByCluster{iClu})
                %                 if S_clu.clusterXPositions(iClu) < posX_thresh, continue; end
                posX_clu1 = S_clu.clusterXPositions(iClu);
                posY_clu1 = S_clu.clusterYPositions(iClu);
                if posX_clu1 < posX_lim(1) || posX_clu1 > posX_lim(2) || posY_clu1 < posY_lim(1) || posY_clu1 > posY_lim(2), continue; end
                plot(ax, vrTime_spk(viSpk1), vrPosY_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5); % - median(vrPosY_spk(viSpk1))
            end %switch
        end
    else
        switch vcMode_drift
            case 'tay', plot3(ax, vrTime_spk, vrAmp_spk, vrPosY_spk, 'o', 'MarkerSize', 5);
            case 'x', plot(ax, vrTime_spk, vrPosX_spk, '.', 'MarkerSize', 5);
            case 'y', plot(ax, vrTime_spk, vrPosY_spk, '.', 'MarkerSize', 5);
            case 'xy', plot(ax, vrPosX_spk, vrPosY_spk, '.', 'MarkerSize', 5);
            case 'xya', plot3(ax, vrPosX_spk, vrPosY_spk, vrAmp_spk, '.', 'MarkerSize', 5);
            case 'z'
            posX_lim = quantile(vrPosX_spk, [.25, .75]);
            posY_lim = quantile(vrPosY_spk, [.35, .65]);
            viSpk_plot = find(vrPosX_spk(:) > posX_lim(1) & vrPosX_spk(:) < posX_lim(2) & vrPosY_spk(:) > posY_lim(1) & vrPosY_spk(:) < posY_lim(2)); % & vrAmp_spk(:) < median(vrAmp_spk));
            plot(ax, vrTime_spk(viSpk_plot), vrPosY_spk(viSpk_plot), '.', 'MarkerSize', 5);
        end %switch
    end

    drawnow;
    title_(ax, sprintf('Spike positions from shank %d (change using "iShank_show" paremter)', iShank_show));
    xlabel(ax, 'Time (s)');
    grid(ax, 'on');
    axis(ax, 'tight');
    switch vcMode_drift
        case 'x', ylabel(ax, 'X position (um)');
        case 'y', ylabel(ax, 'Y position (um)');
        case 'tay'
        ylabel(ax, '1/a (au)');
        zlabel(ax, 'z position (um)');
        case 'xya'
        xlabel(ax, 'x (um)');
        ylabel(ax, 'y (um)');
        zlabel(ax, '1/a (au)');
    end %switch
    mouse_figure(hFig_drift);
    % linkaxes([ax_x, ax_y], 'x');
end %func
