%--------------------------------------------------------------------------
function plot_drift_(P)

    iShank_show = get_set_(P, 'iShank_show', 1);
    vcMode_drift = get_set_(P, 'vcMode_drift', 'z'); % {'tay', 'xy', 'xya', 'x', 'y'}
    vcMode_com = get_set_(P, 'vcMode_com', 'fet'); % {'fet', 'filt', 'raw', 'std'}
    % Compute spike position and plot drift profile
    S0 = load_cached_(P); % load cached data or from file if exists
    [S_clu, spikeSecondarySites, vrTime_spk, spikeSites] = get0_('S_clu', 'spikeSecondarySites', 'spikeTimes', 'spikeSites');
    vrTime_spk = double(vrTime_spk) / P.sRateHz;
    nSites_spk = 1 + P.maxSite*2 - P.nSites_ref;
    miSites_spk = P.miSites(:,spikeSites);

    switch lower(vcMode_com)
        case 'filt'
        tnWav_spk = get_spkwav_(P, 0);
        %mrVp = squeeze_(single(min(tnWav_spk))) .^ 2;
        %         tnWav_spk1 = meanSubt_(single(tnWav_spk),2);
        mrVp = single(squeeze_(max(tnWav_spk) - min(tnWav_spk))) .^ 2;
        case 'filtstd'
        mrVp = squeeze_(var(single(get_spkwav_(P, 0))));
        case 'rawstd', mrVp = squeeze_(var(single(get_spkwav_(P, 1))));
        case 'raw'
        tnWav_raw = get_spkwav_(P, 1);
        mrVp = single(squeeze_(max(tnWav_raw) - min(tnWav_raw))) .^ 2;
        case 'fet'
        spikeFeatures = get_spkfet_(P);
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
    nSpk_thresh_clu = median(S_clu.vnSpk_clu);
    snr_thresh_clu = get_set_(P, 'snr_thresh_clu', quantile(S_clu.vrSnr_clu, .5));
    % posX_thresh = median(S_clu.vrPosX_clu);
    % [vrTime_drift, vrDepth_drift] = drift_track_(S_clu, vrPosY_spk, P);

    if ~isempty(S_clu)
        posX_lim = quantile(S_clu.vrPosX_clu, [.25, .75]);
        posY_lim = quantile(S_clu.vrPosY_clu, [.1, .9]);
        for iClu = 1:S_clu.nClusters
            viSpk1 = S_clu.cviSpk_clu{iClu};
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
                %                 if S_clu.vnSpk_clu(iClu) < nSpk_thresh_clu, continue; end
                if S_clu.vrSnr_clu(iClu) < snr_thresh_clu, continue; end
                %                 if vrA_spk(S_clu.cviSpk_clu{iClu})
                %                 if S_clu.vrPosX_clu(iClu) < posX_thresh, continue; end
                posX_clu1 = S_clu.vrPosX_clu(iClu);
                posY_clu1 = S_clu.vrPosY_clu(iClu);
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
