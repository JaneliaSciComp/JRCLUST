function doPlotFigPos(hFigPos, hClust, hCfg, selected)
    % This also plots cluster position
    [hFigPos, S_fig] = get_fig_cache_('FigPos');

    S_clu1 = hClust.exportUnitInfo(selected(1));
    if numel(selected) > 1
        S_clu2 = hClust.exportUnitInfo(selected(2));
    end

    % plot waveform in space
    if isempty(hFigPos.figData)
        hFigPos.axes();
    else
        hFigPos.cla();
        hFigPos.hold('on');
    end

    plot_unit_(S_clu1, S_fig.hAx, [0 0 0]);
    vrPosXY1 = [hClust.vrPosX_clu(S_clu1.iClu), hClust.vrPosY_clu(S_clu1.iClu)];
    nSpk1 = hClust.vnSpk_clu(S_clu1.iClu);
    if isempty(S_clu2)
        vcTitle = sprintf('Unit %d: %d spikes; (X=%0.1f, Y=%0.1f) [um]', S_clu1.iClu, nSpk1, vrPosXY1);
        try
            vcTitle = sprintf('%s\n%0.1fuVmin, %0.1fuVpp, SNR:%0.1f ISI%%:%2.3f IsoDist:%0.1f L-rat:%0.1f', ...
                vcTitle, S_clu1.uVmin, S_clu1.uVpp, S_clu1.snr, S_clu1.isi_ratio*100, S_clu1.iso_dist,  S_clu1.l_ratio);
        catch
        end
    else
        nSpk2 = hClust.vnSpk_clu(S_clu2.iClu);
        vrPosXY2 = [hClust.vrPosX_clu(S_clu2.iClu), hClust.vrPosY_clu(S_clu2.iClu)] / hCfg.um_per_pix;
        plot_unit_(S_clu2, S_fig.hAx, [1 0 0]);
        vcTitle = sprintf('Unit %d(black)/%d(red); (%d/%d) spikes\n(X=%0.1f/%0.1f, Y=%0.1f/%0.1f) [um]', ...
        S_clu1.iClu, S_clu2.iClu, nSpk1, nSpk2, ...
        [vrPosXY1(1), vrPosXY2(1), vrPosXY1(2), vrPosXY2(2)]);
    end
    title_(S_fig.hAx, vcTitle);
    set(hFigPos, 'UserData', S_fig);
end

%% LOCAL FUNCTIONS
function plot_unit_(S_clu1, hFigPos, vcColor0, hCfg)
    if isempty(S_clu1)
        return;
    end

    [~, S_figWav] = get_fig_cache_('FigWav');
    maxAmp = S_figWav.maxAmp;
    % plot individual unit
    nSamples = size(S_clu1.mrWav_clu, 1);
    vrX = (1:nSamples)'/nSamples;
    vrX([1,end]) = nan; % line break

    if ~isequal(vcColor0, [0 0 0])
        trWav1 = zeros(1,1,0);
    else
        trWav1 = S_clu1.trWav;
    end

    % show example traces
    for iWav = size(trWav1,3):-1:0
        if iWav==0
            mrY1 = S_clu1.mrWav_clu / maxAmp;
            lineWidth=1.5;
            vcColor = vcColor0;
        else
            mrY1 = trWav1(:,:,iWav) / maxAmp;
            lineWidth=.5;
            vcColor = .5*[1,1,1];
        end
        vrX1_site = hCfg.mrSiteXY(S_clu1.viSite, 1) / hCfg.um_per_pix;
        vrY1_site = hCfg.mrSiteXY(S_clu1.viSite, 2) / hCfg.um_per_pix;
        mrY1 = bsxfun(@plus, mrY1, vrY1_site');
        mrX1 = bsxfun(@plus, repmat(vrX, [1, size(mrY1, 2)]), vrX1_site');
        line(mrX1(:), mrY1(:), 'Color', vcColor, 'Parent', hFigPos, 'LineWidth', lineWidth);
    end
    xlabel(hFigPos, 'X pos [pix]');
    ylabel(hFigPos, 'Z pos [pix]');
    grid(hFigPos, 'on');
    xlim_(hFigPos, [min(mrX1(:)), max(mrX1(:))]);
    ylim_(hFigPos, [floor(min(mrY1(:))-1), ceil(max(mrY1(:))+1)]);
end %func
