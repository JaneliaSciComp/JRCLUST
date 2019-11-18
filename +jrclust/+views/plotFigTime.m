function hFigTime = plotFigTime(hFigTime, hClust, hCfg, selected, maxAmp, iSite, channel_idx)
    %DOPLOTFIGTIME Plot features vs. time
    timeLimits = double([0, abs(hClust.spikeTimes(end))/hCfg.sampleRate]);

    % construct plot for the first time
    if ~hFigTime.hasAxes('default')
        hFigTime.addAxes('default');
        hFigTime.axApply('default', @set, 'Position', [0.03 0.2 0.9 0.7], 'XLimMode', 'manual', 'YLimMode', 'manual');

        % trial time indicators
        if isempty(hCfg.trialFile)
            trialTimes=[];
        else
            trialTimes = jrclust.utils.loadTrialFile(hCfg.trialFile);
            if iscell(trialTimes)
                trialTimes = trialTimes{1};
            end
        end
        if ~isempty(trialTimes)
            hFigTime.addPlot('trialTimes',@line,repmat(trialTimes,1,2),[0 abs(maxAmp)],'linewidth',0.1,'color',[0.5 0.7 0.5]);
        elseif ~isempty(hCfg.trialFile)
           warning('Could not load trial times from %s.',hCfg.trialFile);
        end        
        
        % first time
        hFigTime.addPlot('background', @line, nan, nan, 'Marker', '.', 'Color', hCfg.colorMap(1, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.addPlot('foreground', @line, nan, nan, 'Marker', '.', 'Color', hCfg.colorMap(2, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.addPlot('foreground2', @line, nan, nan, 'Marker', '.', 'Color', hCfg.colorMap(3, :), 'MarkerSize', 5, 'LineStyle', 'none');
        hFigTime.axApply('default', @xlabel, 'Time (s)');
        hFigTime.axApply('default', @grid, 'on');

        % rectangle plot
        rectPos = [timeLimits(1), maxAmp, diff(timeLimits), maxAmp];
        hFigTime.addPlot('hRect', @imrect, rectPos);
        hFigTime.plotApply('hRect', @setColor, 'r');
        hFigTime.plotApply('hRect', @setPositionConstraintFcn, makeConstrainToRectFcn('imrect', timeLimits, [-4000 4000]));

        hFigTime.setHideOnDrag('background'); % hide background spikes when dragging
        
        % histogram
        hFigTime.addAxes('histogram');
        hFigTime.axApply('histogram', @set, 'Position', [0.93 0.2 0.06 0.7],'Visible','off'); 
        hist_args = {@stairs, hFigTime.hAxes('histogram'), nan,nan,'LineWidth',1};
        hFigTime.addPlot('background_hist', hist_args{:}, 'color', hCfg.colorMap(1, :));
        hFigTime.addPlot('foreground_hist', hist_args{:}, 'color', hCfg.colorMap(2, :));
        hFigTime.addPlot('foreground_hist2', hist_args{:}, 'color', hCfg.colorMap(3, :));
    end
    
    [bgFeatures, bgTimes] = getFigTimeFeatures(hClust, iSite); % plot background
    [fgFeatures, fgTimes, YLabel] = getFigTimeFeatures(hClust, iSite, selected(1),channel_idx); % plot primary selected cluster

    if numel(selected) == 2
        [fgFeatures2, fgTimes2] = getFigTimeFeatures(hClust, iSite, selected(2));
        figTitle = sprintf('Unit %d (black), Unit %d (red); (press [H] for help)', selected(1), selected(2));
    else
        fgFeatures2 = [];
        fgTimes2 = [];
        figTitle = sprintf('Unit %d (black); (press [H] for help)', selected(1));
    end
   
    binlimits = [min(bgFeatures) max(bgFeatures)]+eps;     
    
    % remove foreground events from background cluster    
    bg_idx = ~ismember(bgTimes,union(fgTimes,fgTimes2));
    bgFeatures = bgFeatures(bg_idx);
    bgTimes = bgTimes(bg_idx);

    vppLim = [0, abs(maxAmp)];
    
    %update scatter plots
    hFigTime.updatePlot('background', bgTimes, bgFeatures);
    hFigTime.updatePlot('foreground', fgTimes, fgFeatures);
    hFigTime.updatePlot('foreground2', fgTimes2, fgFeatures2);
    imrectSetPosition(hFigTime, 'hRect', timeLimits, vppLim);

    % update histograms
    
    n_hist_bins=100; % seems to work nicely
    % add eps so as not to plot background spikes with feature projection of 0. 
    % Sometimes there are a lot of these and they make the other points hard to see in the histogram.
    histcountfun = @(features)histcounts(features,n_hist_bins,'BinLimits',binlimits,'Normalization','probability');
    updateplotfun = @(tag,N,edges)hFigTime.updatePlot(tag,[0 N 0],[edges edges(end)+eps]); % feeding stairs the output of histcounts in this way exactly reproduces the output of matlab histogram, rotated on its side
    [N,edges] = histcountfun(bgFeatures);
    updateplotfun('background_hist',N,edges);
    [N,edges] = histcountfun(fgFeatures);
    updateplotfun('foreground_hist',N,edges);
    [N,edges] = histcountfun(fgFeatures2);    
    updateplotfun('foreground_hist2',N,edges);


%     if isfield(S_fig, 'vhAx_track')
%         toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 0);
%         toggleVisible_({S_fig.hAx, S_fig.hRect, S_fig.hPlot1, S_fig.hPlot2, S_fig.hPlot0}, 1);
%     end
%
    if ~isfield(hFigTime.figData, 'doPlotBG')
        hFigTime.figData.doPlotBG = 1;
    end

    hFigTime.axApply('default', @axis, [timeLimits, vppLim]);
    hFigTime.axApply('default', @title, figTitle);
    hFigTime.axApply('default', @ylabel, YLabel);
end
