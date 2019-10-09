function hFigTime = plotFigTime(hFigTime, hClust, hCfg, selected, maxAmp, iSite, channel_idx)
    persistent linehandle
    %DOPLOTFIGTIME Plot features vs. time
    timeLimits = double([0, abs(hClust.spikeTimes(end))/hCfg.sampleRate]);

    % construct plot for the first time
    if ~hFigTime.hasAxes('default')
        hFigTime.addAxes('default');
        hFigTime.axApply('default', @set, 'Position', [0.03 0.2 0.9 0.7], 'XLimMode', 'manual', 'YLimMode', 'manual');

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
        hist_args = {@histogram, hFigTime.hAxes('histogram'), nan, ...
            'Orientation','horizontal','Normalization','probability',...
            'DisplayStyle','stairs'};
        hFigTime.addPlot('background_hist', hist_args{:}, 'EdgeColor', hCfg.colorMap(1, :));
        hFigTime.addPlot('foreground_hist', hist_args{:}, 'EdgeColor', hCfg.colorMap(2, :));
        hFigTime.addPlot('foreground_hist2', hist_args{:}, 'EdgeColor', hCfg.colorMap(3, :));
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

    binlimits = [min(bgFeatures) max(bgFeatures)];
    
    % remove foreground events from background cluster    
    bg_idx = ~ismember(bgTimes,union(fgTimes,fgTimes2));
    bgFeatures = bgFeatures(bg_idx);
    bgTimes = bgTimes(bg_idx);

    vppLim = [0, abs(maxAmp)];

    hFigTime.updatePlot('background', bgTimes, bgFeatures);
    hFigTime.updatePlot('foreground', fgTimes, fgFeatures);
    hFigTime.updatePlot('foreground2', fgTimes2, fgFeatures2);
    imrectSetPosition(hFigTime, 'hRect', timeLimits, vppLim);

    hFigTime.updateHistogram('background_hist',bgFeatures,[],'BinLimits',binlimits+eps,'NumBins',100);
    hFigTime.updateHistogram('foreground_hist',fgFeatures,[],'BinLimits',binlimits+eps,'NumBins',100);
    hFigTime.updateHistogram('foreground_hist2',fgFeatures2,[],'BinLimits',binlimits+eps,'NumBins',100);


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

    %% add trial time indicators
    if isempty(linehandle)
        trialTimes = loadTrialFile(hCfg.trialFile);
        if ~isempty(trialTimes)
            ax=get(hFigTime.hPlots('foreground'),'parent');
            yl=get(ax,'ylim');
            linehandle = line(ax,repmat(trialTimes{1}(:,1),1,2),yl,'linewidth',0.1,'color',[0 1 0]);
            set(ax,'ylim',yl);
        elseif ~isempty(hCfg.trialFile)
           warning('Could not load trial times from %s.',hCfg.trialFile);
        end
    end
end


function trialTimes = loadTrialFile(trialFile)
    %LOADTRIALFILE Import trial times (in seconds)
    trialTimes = [];

    try
        [~, ~, ext] = fileparts(trialFile);

        if strcmpi(ext, '.mat')
            trialData = load(trialFile);
            fieldNames = fieldnames(trialData);

            trialTimes = trialData.(fieldNames{1});
            if isstruct(trialTimes)
                trialTimes = trialTimes.times;
            end
        elseif strcmpi(ext, '.csv')
            trialTimes = csvread(trialFile);
            if isrow(trialTimes)
                trialTimes = trialTimes(:);
            end
        end
    catch ME
        warning('Could not load trialFile %s: %s', trialFile, ME.message);
    end
end
