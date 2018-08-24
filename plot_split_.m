%--------------------------------------------------------------------------
function [fSplit, vlIn] = plot_split_(S1)
    % find site
    mrPolyPos = getPosition(S1.hPoly);
    site12_show = floor(mean(mrPolyPos));
    site12 = S1.sitesOfInterest(site12_show+1);

    % get amp
    S0 = get(0, 'UserData');
    S_clu = S0.S_clu;
    P = S0.P;
    iClu1 = S0.primarySelectedCluster;
    if ismember(P.displayFeature, {'pca', 'ppca', 'gpca'})
        fWav_raw_show = 0;
    else
        fWav_raw_show = getOr(P, 'fWav_raw_show', 0);
    end
    trWav12 = tnWav2uV_(getSpikeWaveformsSites(S_clu.spikesByCluster{iClu1}, site12, S0, fWav_raw_show), P);
    if diff(site12) == 0, trWav12(:,2,:) = trWav12(:,1,:); end
    vxPoly = (mrPolyPos([1:end,1],1) - site12_show(1)) * S1.maxAmp;
    vyPoly = (mrPolyPos([1:end,1],2) - site12_show(2)) * S1.maxAmp;
    switch lower(P.displayFeature)
        case {'vpp', 'vmin', 'vmax'}
        mrAmin12 = abs(squeeze_(min(trWav12)))';
        mrAmax12 = abs(squeeze_(max(trWav12)))';
        vyPlot = mrAmin12(:,2);
        vcYlabel = sprintf('Site %d (min amp)', site12(2));
        if site12(2) > site12(1)
            vxPlot = mrAmin12(:,1);
            vcXlabel = sprintf('Site %d (min amp)', site12(1));
        else
            vxPlot = mrAmax12(:,1);
            vcXlabel = sprintf('Site %d (max amp)', site12(1));
            vxPoly = vxPoly; %max amp are scaled half
        end

        case {'cov', 'spacetime'}
        [mrAmin12, mrAmax12] = calc_cov_spk_(S_clu.spikesByCluster{iClu1}, site12);
        [mrAmin12, mrAmax12] = multifun_(@(x)abs(x'), mrAmin12, mrAmax12);
        vyPlot = mrAmin12(:,2);
        vcYlabel = sprintf('Site %d (cov1)', site12(2));
        if site12(2) > site12(1)
            vxPlot = mrAmin12(:,1);
            vcXlabel = sprintf('Site %d (cov1)', site12(1));
        else
            vxPlot = mrAmax12(:,1);
            vcXlabel = sprintf('Site %d (cov2)', site12(1));
        end

        case {'pca', 'ppca', 'gpca'}
        if strcmpi(P.displayFeature, 'ppca')
            [mrPv1, mrPv2] = pca_pv_clu_(site12, S0.primarySelectedCluster, S0.secondarySelectedCluster);
            [mrAmin12, mrAmax12] = pca_pc_spk_(S_clu.spikesByCluster{iClu1}, site12, mrPv1, mrPv2);
        else
            [mrAmin12, mrAmax12] = pca_pc_spk_(S_clu.spikesByCluster{iClu1}, site12);
        end
        [mrAmin12, mrAmax12] = multifun_(@(x)abs(x'), mrAmin12, mrAmax12);
        vyPlot = mrAmin12(:,2);
        vcYlabel = sprintf('Site %d (PC1)', site12(2));
        if site12(2) > site12(1)
            vxPlot = mrAmin12(:,1);
            vcXlabel = sprintf('Site %d (PC1)', site12(1));
        else
            vxPlot = mrAmax12(:,1);
            vcXlabel = sprintf('Site %d (PC2)', site12(1));
        end

        otherwise
        error('plot_split: featureShow: not implemented');
        %         vxPoly = (mrPolyPos([1:end,1],1) - site12_show(1)) * S1.maxAmp;
        %         vyPoly = (mrPolyPos([1:end,1],2) - site12_show(2)) * S1.maxAmp;
        %         trFet12 = trFet_([], site12, S_clu.spikesByCluster{iClu1});
        %         vyPlot = squeeze_(trFet12(1, 2, :));
        %         vcYlabel = sprintf('Site %d (%s1)', site12(2), P.feature);
        %         if site12(2) > site12(1)
        %             vxPlot = squeeze_(trFet12(1, 1, :));
        %             vcXlabel = sprintf('Site %d (%s1)', site12(1), P.feature);
        %         else
        %             vxPlot = squeeze_(trFet12(min(2,P.nPcPerChan), 1, :));
        %             vcXlabel = sprintf('Site %d (%s2)', site12(1), P.feature);
        %         end
    end

    vlIn = inpolygon(vxPlot, vyPlot, vxPoly, vyPoly);

    % Plot temporary figure (auto-close)
    hFig = figure(10221); clf;
    resizeFigure(hFig, [0 0 .5 1]);
    subplot(2,2,[1,3]); hold on;
    line(vxPlot, vyPlot, 'Color', P.mrColor_proj(2,:), 'Marker', 'o', 'MarkerSize', 2, 'LineStyle', 'none');
    hPlot = line(vxPlot(vlIn), vyPlot(vlIn), 'Color', P.mrColor_proj(3,:), 'Marker', 'o', 'MarkerSize', 2, 'LineStyle', 'none');
    % plot(vxPoly, vyPoly, 'b+-'); %boundary
    title(sprintf('Cluster %d (%d spikes)', iClu1, S_clu.nSpikesPerCluster(iClu1)));
    xlabel(sprintf('Site %d', site12(1)));
    ylabel(vcYlabel);   xlabel(vcXlabel);
    grid on;

    % user edit polygon
    hPoly = impoly_(gca, [vxPoly(:), vyPoly(:)]);
    hFunc = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(hPoly, hFunc);

    hMsgbox = msgbox_('Press OK after adjusting polygon', 0);
    uiwait(hMsgbox);

    vlIn = poly_mask_(hPoly, vxPlot, vyPlot);
    set(hPlot, 'XData', vxPlot(vlIn), 'YData',vyPlot(vlIn));
    % if P.fWav_raw_show
    %     trWav12_raw = tnWav2uV_(getSpikeWaveformsSites(S_clu.spikesByCluster{iClu1}, site12, 1));
    %     mrWavX = squeeze_(trWav12_raw(:, 1, :));
    %     mrWavY = squeeze_(trWav12_raw(:, 2, :));
    % else
    mrWavX = squeeze_(trWav12(:, 1, :));
    mrWavY = squeeze_(trWav12(:, 2, :));
    % end
    vrT = (P.spkLim(1):P.spkLim(end)) / P.sampleRateHz * 1000;
    viIn = randomSubsample(find(vlIn), P.nSpk_show);
    viOut = randomSubsample(find(~vlIn), P.nSpk_show);

    if isempty(viIn) || isempty(viOut)
        fSplit = 0;
        tryClose(hFig);
        return;
    end

    subplot 222; hold on;
    plot(vrT, mrWavX(:,viOut), 'k', vrT, mrWavX(:,viIn), 'r');
    title(vcXlabel); ylabel('Voltage (\muV)'); xlabel('Time (ms)');
    grid on;

    subplot 224; hold on;
    plot(vrT, mrWavY(:,viOut), 'k', vrT, mrWavY(:,viIn), 'r');
    title(vcYlabel); ylabel('Voltage (\muV)'); xlabel('Time (ms)');
    grid on;

    if strcmpi(userDialog('Split?', 'Confirmation', 'Yes', 'No', 'Yes'), 'Yes')
        fSplit = 1;
    else
        fSplit = 0;
    end
    tryClose(hFig);
end %func
