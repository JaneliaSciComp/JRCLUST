%--------------------------------------------------------------------------
function [fSplit, vlIn] = plot_split_(S1)
    % find site
    mrPolyPos = getPosition(S1.hPoly);
    site12_show = floor(mean(mrPolyPos));
    site12 = S1.viSites_show(site12_show+1);

    % get amp
    S0 = get(0, 'UserData');
    hClust = S0.hClust;
    P = S0.P;
    iClu1 = S0.iCluCopy;
    if ismember(P.vcFet_show, {'pca', 'ppca', 'gpca'})
        showRaw = 0;
    else
        showRaw = get_set_(P, 'showRaw', 0);
    end
    trWav12 = jrclust.utils.filtTouV(jrclust.utils.getSampledWindows(hClust.spikesByCluster{iClu1}, site12, S0, showRaw), P);
    if diff(site12) == 0, trWav12(:,2,:) = trWav12(:,1,:); end
    vxPoly = (mrPolyPos([1:end,1],1) - site12_show(1)) * S1.maxAmp;
    vyPoly = (mrPolyPos([1:end,1],2) - site12_show(2)) * S1.maxAmp;
    switch lower(P.vcFet_show)
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
        [mrAmin12, mrAmax12] = getSpikeCov(hClust, hClust.spikesByCluster{iClu1}, site12);
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
        if strcmpi(P.vcFet_show, 'ppca')
            %[prVecs1, prVecs2] = pca_pv_clu_(site12, S0.iCluCopy, S0.iCluPaste);
            [prVecs1, prVecs2] = jrclust.features.getPVClusters(hClust, site12, S0.iCluCopy, S0.iCluPaste);
            [mrAmin12, mrAmax12] = jrclust.features.pcProjectSpikes(hClust.spikesByCluster{iClu1}, site12, prVecs1, prVecs2);
        else
            sampledWindows = permute(jrclust.utils.getSampledWindows(hClust, hClust.spikesByCluster{iClu1}, site12, 0), [1, 3, 2]); % nSamples x nSpikes x nSites
            [prVecs1, prVecs2, prVecs3] = jrclust.features.getPVSpikes(sampledWindows);
            [mrAmin12, mrAmax12] = jrclust.features.pcProjectSpikes(sampledWindows, prVecs1, prVecs2, prVecs3);
            %[mrAmin12, mrAmax12] = jrclust.features.pcProjectSpikes(hClust.spikesByCluster{iClu1}, site12);
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
        error('plot_split: vcFetShow: not implemented');
        %         vxPoly = (mrPolyPos([1:end,1],1) - site12_show(1)) * S1.maxAmp;
        %         vyPoly = (mrPolyPos([1:end,1],2) - site12_show(2)) * S1.maxAmp;
        %         trFet12 = trFet_([], site12, hClust.spikesByCluster{iClu1});
        %         vyPlot = squeeze_(trFet12(1, 2, :));
        %         vcYlabel = sprintf('Site %d (%s1)', site12(2), P.vcFet);
        %         if site12(2) > site12(1)
        %             vxPlot = squeeze_(trFet12(1, 1, :));
        %             vcXlabel = sprintf('Site %d (%s1)', site12(1), P.vcFet);
        %         else
        %             vxPlot = squeeze_(trFet12(min(2,P.nPcPerChan), 1, :));
        %             vcXlabel = sprintf('Site %d (%s2)', site12(1), P.vcFet);
        %         end
    end

    vlIn = inpolygon(vxPlot, vyPlot, vxPoly, vyPoly);

    % Plot temporary figure (auto-close)
    hFig = figure(10221); clf;
    resize_figure_(hFig, [0 0 .5 1]);
    subplot(2,2,[1,3]); hold on;
    line(vxPlot, vyPlot, 'Color', P.mrColor_proj(2,:), 'Marker', 'o', 'MarkerSize', 2, 'LineStyle', 'none');
    hPlot = line(vxPlot(vlIn), vyPlot(vlIn), 'Color', P.mrColor_proj(3,:), 'Marker', 'o', 'MarkerSize', 2, 'LineStyle', 'none');
    % plot(vxPoly, vyPoly, 'b+-'); %boundary
    title(sprintf('Cluster %d (%d spikes)', iClu1, hClust.vnSpk_clu(iClu1)));
    xlabel(sprintf('Site %d', site12(1)));
    ylabel(vcYlabel);   xlabel(vcXlabel);
    grid on;

    % user edit polygon
    hPoly = impoly_(gca, [vxPoly(:), vyPoly(:)]);
    hFunc = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(hPoly, hFunc);

    hMsgbox = jrclust.utils.qMsgBox('Press OK after adjusting polygon', 0);
    uiwait(hMsgbox);

    vlIn = poly_mask_(hPoly, vxPlot, vyPlot);
    set(hPlot, 'XData', vxPlot(vlIn), 'YData',vyPlot(vlIn));
    % if P.showRaw
    %     trWav12_raw = jrclust.utils.filtTouV(jrclust.utils.getSampledWindows(hClust.spikesByCluster{iClu1}, site12, 1));
    %     mrWavX = squeeze_(trWav12_raw(:, 1, :));
    %     mrWavY = squeeze_(trWav12_raw(:, 2, :));
    % else
    mrWavX = squeeze_(trWav12(:, 1, :));
    mrWavY = squeeze_(trWav12(:, 2, :));
    % end
    vrT = (P.spkLim(1):P.spkLim(end)) / P.sRateHz * 1000;
    viIn = randomSelect_(find(vlIn), P.nSpk_show);
    viOut = randomSelect_(find(~vlIn), P.nSpk_show);

    if isempty(viIn) || isempty(viOut)
        fSplit = 0;
        jrclust.utils.tryClose(hFig);
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

    if strcmpi(questdlg_('Split?', 'Confirmation', 'Yes', 'No', 'Yes'), 'Yes')
        fSplit = 1;
    else
        fSplit = 0;
    end
    jrclust.utils.tryClose(hFig);
end %func
